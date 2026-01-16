% Subject folders
subjects = {'EC304','EC288','PR05','PR06','BJH058','DP01'};

% Save starting directory
baseDir = pwd;

% Initialize storage
results1 = [];
results2 = [];

alpha_reward_store = [];
alpha_punish_store = [];
lambda_store = [];
mu_store = [];

for n = 1:numel(subjects)

    % Change into subject directory
    subjDir = fullfile(baseDir, subjects{n});
    cd(subjDir);

    % Load subject-specific behavior
    load('behavior.mat');

    % Run model
    [rho_conflict, pval_conflict, rho_rt, pval_rt, ...
     alpha_reward, alpha_punish, lambda, mu] = ...
        plotSubjectData(n, NaN, NaN, ...
                        reward_trial, punish_trial, ...
                        decision, conflict_trial_type, ...
                        reaction_time, reward_trial_type, ...
                        punishment_trial_type, p_approach_trial_type, ...
                        conflict_trial);


    % Store results
    results1 = [results1; n, rho_conflict, pval_conflict];
    results2 = [results2; n, rho_rt, pval_rt];

    alpha_reward_store = [alpha_reward_store alpha_reward];
    alpha_punish_store = [alpha_punish_store alpha_punish];
    lambda_store      = [lambda_store lambda];

    % Return to base directory
    cd(baseDir);

end

% Display results
disp('Conflict vs |Value| (row = subject index)');
disp(results1);

disp('RT vs |Value| (row = subject index)');
disp(results2);
%% Extended Data Fig 1f  parameter fits

a_r = alpha_reward_store;
a_p   = alpha_punish_store;
lamb = lambda_store;

labels = {'a_r','a_p','lambda'};
x_ar = 1 * ones(size(a_r));
x_ap   = 2 * ones(size(a_p));
x_lamb = 3 * ones(size(lamb));

% Shared beeswarm so points don't overlap within each column
[x_ar_sw, y_ar_sw] = beeswarm1D(x_ar, a_r, 0.35);
[x_ap_sw,   y_ap_sw]   = beeswarm1D(x_ap,   a_p,   0.35);
[x_lamb_sw,   y_lamb_sw]   = beeswarm1D(x_lamb,   lamb,   0.35);

figure; hold on

% Predefine colors
gray_levels = (0:5) * 0.15;           % one per subject

for n = 1:6
    g = gray_levels(n);

    % Alpha-reward (subject-specific grayscale)
    scatter(x_ar_sw(n), y_ar_sw(n), 70, ...
        [g g g], 'filled');

    % Alpha-punish
    scatter(x_ap_sw(n), y_ap_sw(n), 70, ...
        [g g g], 'filled');

    % Lambda
    scatter(x_lamb_sw(n), y_lamb_sw(n), 70, ...
        [g g g], 'filled');
end


box off
set(gca,'TickDir','out')

xlim([0.5 3.5])
xticks([1 2 3])
ylim([0 3])
yline([1 1])
xticklabels(labels)
ylabel('Parameter fit (a.u.)')
yline(0,'k-');
box off; set(gca,'TickDir','out');
title('Subject-level parameter fits')
hold off
%% Summary Plot Fig 1

% results1(:,2) = rho_conflict, results2(:,2) = rho_rt
rho_conf = results1(:,2);
rho_rt   = results2(:,2);

labels = {'Conflict \rho','RT \rho'};
x_conf = 1 * ones(size(rho_conf));
x_rt   = 2 * ones(size(rho_rt));

% Shared beeswarm so points don't overlap within each column
[x_conf_sw, y_conf_sw] = beeswarm1D(x_conf, rho_conf, 0.18);
[x_rt_sw,   y_rt_sw]   = beeswarm1D(x_rt,   rho_rt,   0.18);

figure; hold on
scatter(x_conf_sw, y_conf_sw, 70, 'filled');  % conflict points
scatter(x_rt_sw,   y_rt_sw,   70, 'filled');  % RT points


xlim([0.5 2.5])
xticks([1 2])
ylim([-1 1])
xticklabels(labels)
ylabel('\rho (Spearman)')
yline(0,'k-');
box off; set(gca,'TickDir','out');
title('Subject-level Spearman \rho: Conflict vs RT')
hold off

%% Adapted plotSubjectData function using prospect theory model
function [rho_conflict, pval_conflict, rho_rt, pval_rt,alpha_reward_opt,alpha_punish_opt,lambda_opt,mu_opt] = ...
    plotSubjectData(subjectID, reward_coeff, punish_coeff, rewards_trials, punishments_trials, ...
                    decision, conflict_all, reaction_time_all, reward_trial_types, punishment_trial_types, p_approach,conflict_trial)
    filter = find(decision<2);
    % Use the passed-in trial data for plotting
    rewards = rewards_trials(filter);
    punishments = punishments_trials(filter);
    decisions = decision(filter);
    %conflict = conflict_all(filter);
    reaction_times = reaction_time_all(filter);
    %conflict_trial = conflict_trial(filter);
    
    % For grid-based estimation, use the unique reward and punishment sizes:
    reward_sizes = 0:1:7;
    punishment_sizes = 0:1:4;
    
    % Build a matrix (grid) for the observed p_approach values.
    % Rows correspond to punishment_sizes and columns to reward_sizes.
    mat = NaN(length(punishment_sizes), length(reward_sizes));
    for i = 1:length(punishment_sizes)
        for j = 1:length(reward_sizes)
            idx = find(punishment_trial_types == punishment_sizes(i) & reward_trial_types == reward_sizes(j));
            if ~isempty(idx)
                mat(i, j) = p_approach(idx);  % If more than one index exists, this uses the available data.
            end
        end
    end

    % Flatten the grid and remove NaNs (to build model fitting data)
    [X_grid, Y_grid] = meshgrid(reward_sizes, punishment_sizes);
    X_flat = X_grid(:);
    Y_flat = Y_grid(:);
    P_flat = mat(:);
    valid_idx = ~isnan(P_flat);
    X_flat = X_flat(valid_idx);
    Y_flat = Y_flat(valid_idx);
    P_flat = P_flat(valid_idx);

    %%% Prospect Theory Model Fitting Section
    % Known outcome probabilities (as provided):
    p_reward = 0.42;      % (0.6 * 0.7)
    p_punishment = 0.4;
    
    % Estimate parameters for the model using the grid data.
    % The cost function (costfun) computes the negative log-likelihood:
    %   params = [mu, alpha_reward, alpha_punish, lambda]
    objfun = @(params) costfun(params, X_flat, Y_flat, P_flat, p_reward, p_punishment);
    init_params = [1, 1, 1, 1];  % initial guesses for [mu, α_reward, α_punish, λ]
    options = optimset('Display','off','TolFun',1e-8);
    [opt_params, ~] = fminsearch(objfun, init_params, options);
    
    % Unpack optimized parameters:
    mu_opt = opt_params(1);
    alpha_reward_opt = opt_params(2);
    alpha_punish_opt = opt_params(3);
    lambda_opt = opt_params(4);
    
    % Compute the decision value (V) on the grid:
    V_grid = p_reward * (X_grid.^alpha_reward_opt) - p_punishment * lambda_opt * (Y_grid.^alpha_punish_opt);
    
    %%% Plotting Section
    % SUBPLOT 2: Plot the decision function (image) and overlay the decision boundary (V = 0)
    subplot(1, 6, 2)
    imagesc(reward_sizes, punishment_sizes, V_grid)
    colormap bone
    colorbar
    hold on;
    % The decision boundary is where V = 0 (p_pred = 0.5)
    contour(reward_sizes, punishment_sizes, V_grid, [0 0], 'LineColor', 'k', 'LineWidth', 2);
    set(gca, 'YDir', 'normal');
    xlabel('Reward size');
    ylabel('Punishment size');
    title('Decision boundary (V=0)');
    box off; set(gca, 'TickDir', 'out');
    hold off;

    % SUBPLOT 1: TRUE 2D beeswarm (shared across decisions → no overlap)
    subplot(1, 6, 1)

    % Build a single swarm for ALL trials
    [x_sw, y_sw] = beeswarm2D_on_grid(rewards, punishments, 0.5);

    % Plot by decision type (same coordinates, different markers)
    scatter(x_sw(decisions==1), y_sw(decisions==1), ...
        50, 'o', 'MarkerEdgeColor', 'b'); hold on

    scatter(x_sw(decisions==0), y_sw(decisions==0), ...
        50, '+', 'MarkerEdgeColor', 'r');

    % Overlay decision boundary
    contour(reward_sizes, punishment_sizes, V_grid, [0 0], ...
        'LineColor', 'k', 'LineWidth', 2);

    xlim([-0.5 7.5]); ylim([-0.5 4.5]);
    xlabel('Reward size'); ylabel('Punishment size');
    title('Trial data & decision boundary (2D beeswarm)');
    box off; set(gca, 'TickDir', 'out');
    hold off


    
    % SUBPLOT 3: Conflict vs. |Value|
    subplot(1, 6, 3)
    % Compute trial-level decision value using the prospect theory model
    % (Using individual trial offers, with the same estimated exponents)
    Value = p_reward * (reward_trial_types.^alpha_reward_opt) - p_punishment * lambda_opt * (punishment_trial_types.^alpha_punish_opt);
    %decision_function_conflict = abs(Val_apr_conflict);
    % Plot the data points
    plot(abs(Value), conflict_all, 'k.', 'MarkerSize', 20);
    hold on;
    % Fit a linear regression model and plot the best-fit line
    [rho_conflict, pval_conflict] = corr(abs(Value)', conflict_all','Type','Spearman');
    mdl_conflict = fitlm(abs(Value), conflict_all);
    beta_conflict = mdl_conflict.Coefficients.Estimate(2);
    x_fit = linspace(min(abs(Value)), max(abs(Value)), 100);
    y_fit = mdl_conflict.Coefficients.Estimate(1) + beta_conflict*x_fit;
    plot(x_fit, y_fit, 'r-', 'LineWidth', 2);
    text(max(x_fit), max(y_fit), sprintf('p = %.3f, rho = %.3f', pval_conflict, rho_conflict), ...
         'FontSize', 12, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right');
    xlabel('|Decision Value|'); ylabel('Conflict');
    title('Conflict vs. |Value|');
    box off; set(gca, 'TickDir', 'out');
    hold off;

    % SUBPLOT 4: Reaction Time vs. |Value|
    subplot(1, 6, 4)
    % Compute trial-level decision value for reaction time analysis
    Val_apr_rt = p_reward * (rewards.^alpha_reward_opt) - p_punishment * lambda_opt * (punishments.^alpha_punish_opt);
    decision_function_rt = abs(Val_apr_rt);
    scatter(decision_function_rt + 0.01*randn(length(decision_function_rt),1), reaction_times, 5, 'k.');
    hold on;
    [rho_rt, pval_rt] = corr(decision_function_rt(decisions < 2), reaction_times(decisions < 2),'Type','Spearman');
    mdl_rt = fitlm(decision_function_rt, reaction_times);
    beta_rt = mdl_rt.Coefficients.Estimate(2);
    x_fit_rt = linspace(min(decision_function_rt), max(decision_function_rt), 100);
    y_fit_rt = mdl_rt.Coefficients.Estimate(1) + beta_rt*x_fit_rt;
    plot(x_fit_rt, y_fit_rt, 'r-', 'LineWidth', 2);
    text(max(x_fit_rt), max(y_fit_rt), sprintf('p = %.3f, rho = %.3f', pval_rt, rho_rt), ...
         'FontSize', 12, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right');
    xlabel('|Decision Value|'); ylabel('Reaction time (s)');
    title('RT vs. |Value|');
    box off; set(gca, 'TickDir', 'out');
    hold off;

    % SUBPLOT 5: Conflict vs. Value
    subplot(1, 6, 5)
    % Compute trial-level decision value using the prospect theory model
    % (Using individual trial offers, with the same estimated exponents)
    Value = p_reward * (reward_trial_types.^alpha_reward_opt) - p_punishment * lambda_opt * (punishment_trial_types.^alpha_punish_opt);
    %decision_function_conflict = abs(Val_apr_conflict);
    % Plot the data points
    plot(Value, conflict_all, 'k.', 'MarkerSize', 20);
    lim = max(abs([min(Value) max(Value)]));
    xlim([-lim-0.5 lim+0.5])
    hold on;
    title('Conflict vs. Value');
    box off; set(gca, 'TickDir', 'out');
    hold off;

    % SUBPLOT 6: Reaction Time vs. Value
    subplot(1, 6, 6)
    % Compute trial-level decision value for reaction time analysis
    Val_apr_rt = p_reward * (rewards.^alpha_reward_opt) - p_punishment * lambda_opt * (punishments.^alpha_punish_opt);
    decision_function_rt = Val_apr_rt;
    scatter(decision_function_rt + 0.01*randn(length(decision_function_rt),1), reaction_times, 5, 'k.');
    hold on;
    xlabel('|Decision Value|'); ylabel('Reaction time (s)');
    xlim([-lim-0.5 lim+0.5])
    title('RT vs. Value');
    box off; set(gca, 'TickDir', 'out');
    hold off;

end

%% Cost Function for Prospect Theory Model
function err = costfun(params, X, Y, Pobs, p_reward, p_punishment)
    % params = [mu, alpha_reward, alpha_punish, lambda]
    mu = params(1);
    alpha_reward = params(2);
    alpha_punish = params(3);
    lambda = params(4);
    
    % Compute the approach value for each observation:
    V = p_reward * (X.^alpha_reward) - p_punishment * lambda * (Y.^alpha_punish);
    
    % Predicted probability via logistic function:
    p_pred = 1./(1+exp(-mu*V));
    
    % Clip probabilities to avoid log(0):
    epsilon = 1e-10;
    p_pred = min(max(p_pred, epsilon), 1-epsilon);
    
    % Negative log-likelihood:
    LL = sum(Pobs.*log(p_pred) + (1-Pobs).*log(1-p_pred));
    err = -LL;  % Minimization target
end

function [x_sw, y_sw] = beeswarm2D_on_grid(x, y, radius)
%BEESWARM2D_ON_GRID
%   Shared 2D beeswarm within each discrete (x,y) grid cell.
%   All points in a cell are packed together so different markers never overlap.
%
%   x,y    : discrete coordinates (e.g. reward & punishment sizes)
%   radius : max offset from cell center (≈ 0.15–0.25)

    x = x(:); y = y(:);
    x_sw = x; y_sw = y;

    % Unique grid cells
    cells = unique([x y], 'rows');

    % Spiral packing parameters
    step   = radius / 3;
    golden = pi * (3 - sqrt(5));  % golden angle

    for c = 1:size(cells,1)
        cx = cells(c,1);
        cy = cells(c,2);

        idx = find(x==cx & y==cy);
        n = numel(idx);
        if n <= 1, continue; end

        % Stable ordering (deterministic)
        idx = sort(idx);

        for k = 1:n
            i = k - 1;  % i=0,1,2,...
            if i == 0
                dx = 0; dy = 0;
            else
                r = step * sqrt(i);
                theta = i * golden;
                dx = r * cos(theta);
                dy = r * sin(theta);

                % Clip to radius
                rr = hypot(dx, dy);
                if rr > radius
                    s = radius / rr;
                    dx = dx * s;
                    dy = dy * s;
                end
            end

            x_sw(idx(k)) = cx + dx;
            y_sw(idx(k)) = cy + dy;
        end
    end
end
function [x_sw, y_sw] = beeswarm1D(x_cat, y, width)
%BEESWARM1D Beeswarm-style packing along x for a categorical swarm plot.
%   x_cat : constant category positions (e.g., all 1s or all 2s)
%   y     : values to plot on y axis
%   width : max horizontal spread around category center

    x_cat = x_cat(:);
    y     = y(:);

    x_sw = x_cat;
    y_sw = y;

    cats = unique(x_cat)';
    step = width/4;                 % horizontal spacing
    max_k = floor(width/step);

    for c = cats
        idx = find(x_cat == c);
        if numel(idx) <= 1, continue; end

        % Sort by y so swarm looks neat
        [~, ord] = sort(y(idx), 'ascend');
        idx = idx(ord);

        n = numel(idx);

        % symmetric offsets: 0, +1, -1, +2, -2, ...
        k = zeros(n,1);
        for i = 1:n
            if i == 1
                k(i) = 0;
            else
                t = ceil((i-1)/2);
                k(i) = t * (1 - 2*mod(i,2)); % even->+t, odd->-t
            end
        end

        k = max(min(k, max_k), -max_k);
        x_sw(idx) = c + step*k;
    end
end


