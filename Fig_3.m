
%% First organize electrodes that are all on on hemisphere of each subject's brain
% 1 = right side; 2 = left side
% SUBJECT 1
subject(1).array=zeros(length(subject(1).electrode),1);
for i=3:7
subject(1).array(i) = 1;
end
for i=10:16
subject(1).array(i) = 2;
end

% SUBJECT 2
subject(2).array=zeros(length(subject(2).electrode),1);
for i=4:9
subject(2).array(i) = 1;
end
for i=10:17
subject(2).array(i) = 2;
end

% SUBJECT 3
subject(3).array=zeros(length(subject(3).electrode),1);
for i=4:15
subject(3).array(i) = 1;
end
for i=16:24 % this is the only subject who had 2 arrays on one side hence this one also getting labeled "1"
subject(3).array(i) = 1;
end
for i=28:37
subject(3).array(i) = 2;
end

% SUBJECT 4
subject(4).array=zeros(length(subject(4).electrode),1);
for i=4:14
subject(4).array(i) = 1;
end
for i=18:22
subject(4).array(i) = 2;
end

% SUBJECT 5
subject(5).array=zeros(length(subject(5).electrode),1);
for i=9:17
subject(5).array(i) = 1;
end

% SUBJECT 6
subject(6).array=zeros(length(subject(6).electrode),1);
for i=3:10
subject(6).array(i) = 1;
end
for i=13:20
subject(6).array(i) = 2;
end


%% noise correlations: normalized by value and decision
 

bound1 = 5700; bound2 = 6000;
total_comparisons = 0;
neg_corrs=[];noise_correlation_store = [];null_corrs=[];


    neg_corrs=[];noise_correlation_store = [];all_pairs = [];  % will store ALL tested electrode pairs with p and corr

    for array = 1:2
        % Initialize noise_electrode cell array
        for n = 1:6
            clear I;
            clear coordinate_med;clear coordinate_trans;clear noise_electrode;
            electrodes = find(subject(n).array==array);
            clear noise_electrode;


            % Compute noise correlation matrix and p-values
            num_electrodes = length(electrodes);

            counter = 1;

            for i = 1:num_electrodes
                e = electrodes(i);

                % --- Use all decision 0/1 trials for this subject ---
                which_trials_all = find(subject(n).decision <2 & subject(n).electrode(e).trigger(3).good_trials==1 & subject(n).electrode(e).trigger(1).good_trials==1);
                if isempty(which_trials_all)
                    continue;
                end

                % --- Get value per trial via trial type ---
                % value_trial_type: value for each trial type (1..X)
                % trial_type_trial: trial type index (1..X) for each trial (1..N)
                trial_type_trial = subject(n).trial_type_trial;       % [Ntrials x 1]
                value_trial_type = subject(n).value_trial_type;       % [X x 1]

                val_per_trial = value_trial_type(trial_type_trial);   % [Ntrials x 1]

                % Restrict to decision 0/1 trials
                dec_all = subject(n).decision(which_trials_all);      % decisions for selected trials
                val_all = val_per_trial(which_trials_all);           % values for selected trials

                % --- Baseline-normalized HG in [bound1:bound2] ---
                HG_dec = subject(n).electrode(e).trigger(3).high_gamma_mat;  % [trials x time]
                BL     = subject(n).electrode(e).trigger(1).high_gamma_mat;  % [trials x time]

                % Baseline over 5500:6000
                baseline_all = mean(BL(which_trials_all, 5500:6000), 2);     % [nSelTrials x 1]

                % Signal in analysis window
                sig_all = HG_dec(which_trials_all, bound1:bound2);           % [nSelTrials x T]

                % If these are gpuArrays, bring them back to CPU
                if isa(baseline_all, 'gpuArray')
                    baseline_all = gather(baseline_all);
                end
                if isa(sig_all, 'gpuArray')
                    sig_all = gather(sig_all);
                end

                % Normalize by baseline: (HG - BL) ./ BL
                % (implicit expansion OK in recent MATLAB; use bsxfun if older)
                norm_all = (sig_all - baseline_all) ./ baseline_all;         % [nSelTrials x T]

                % --- Condition labels: (decision, value) per trial ---
                cond_mat = [dec_all(:), val_all(:)];                        % [nSelTrials x 2]
                [~, ~, cond_id] = unique(cond_mat, 'rows');                 % cond_id: 1..nConds

                nConds    = max(cond_id);
                noise_all = zeros(size(norm_all));                          % [nSelTrials x T]

                % Subtract mean trace for each (decision, value) condition
                for c = 1:nConds
                    idx_c = (cond_id == c);                  % trials in this (decision,value) combo
                    mu_c  = mean(norm_all(idx_c, :), 1);     % [1 x T] condition-mean trace
                    noise_all(idx_c, :) = norm_all(idx_c, :) - mu_c;
                end

% --- Normality test on per-trial mean noise (scalar per trial) ---
noise_trial_mean = mean(noise_all, 2);     % [nSelTrials x 1]

% Lilliefors test for normality
[h_lillie, p_lillie] = lillietest(noise_trial_mean);

% store per electrode if you want
lillie_h(counter) = h_lillie;
lillie_p(counter) = p_lillie;


                % Store this electrode's noise (trials x time)
                noise_electrode{counter} = noise_all;

                % Coordinates unchanged
                coordinate_med(counter)   = subject(n).electrode(e).med_coor;
                coordinate_trans(counter) = subject(n).electrode(e).trans_coor;

                counter = counter + 1;
            end


            noise_correlation_matrix = zeros(counter-1, counter-1);
            pval_matrix = zeros(counter-1, counter-1);

            if exist('coordinate_med')
                [B,I] = sort(coordinate_med);

                for i = 1:counter-1
                    for j = i:counter-1
                        if i ==j | abs(i - j)<1
                            noise_correlation_matrix(I(i), I(j)) = 1;
                            pval_matrix(I(i), I(j)) = 1;
                        else

                            % Compute Pearson correlation and p-value
                            [corr_value, p_value] = calculate_corr(noise_electrode{I(i)}, noise_electrode{I(j)});

                            % Store in matrices
                            noise_correlation_matrix(i,j) = corr_value;
                            noise_correlation_matrix(j,i) = corr_value; % Symmetric matrix
                            if corr_value < 0
                                neg_corrs = [neg_corrs; n electrodes(I(i)) electrodes(I(j)) B(i) B(j) coordinate_trans(I(i)) coordinate_trans(I(j)) p_value corr_value];
                            else

                            end

                            null_corrs = [null_corrs; n electrodes(I(i)) electrodes(I(j)) B(i) B(j) coordinate_trans(I(i)) coordinate_trans(I(j)) corr_value];
                            all_pairs = [all_pairs; n array electrodes(I(i)) electrodes(I(j)) B(i) B(j) coordinate_trans(I(i)) coordinate_trans(I(j)) p_value corr_value];

                            total_comparisons = total_comparisons + 1;

                            pval_matrix(i,j) = p_value;
                            pval_matrix(j,i) = p_value; % Symmetric matrix
                        end
                    end
                end

                % Convert the matrix into a single column vector
                vector_column1 = noise_correlation_matrix(:);

                % Create a second column with the value of n repeated
                vector_column2 = repmat(n, length(vector_column1), 1);

                % Combine into a two-column matrix
                result_matrix = [vector_column1, vector_column2];

                noise_correlation_store = [noise_correlation_store;result_matrix];

                noise_correlation_matrices{n} = noise_correlation_matrix;
                pval_matrices{n} = log10(pval_matrix);

                % Subplot 1: Noise correlation matrix
                %if n ==4 && array == 2 %plotting example matrix
                subplot(6,3,(n-1)*3+array)
                imagesc(noise_correlation_matrices{n}); % Display correlation matrix
                colorbar;
                colormap gray
                caxis([-0.05 0])
                all_electrodes = electrodes(I);
                num_electrodes = length(all_electrodes);
                set(gca, 'XTick', 1:num_electrodes, 'XTickLabel', all_electrodes);
                set(gca, 'YTick', 1:num_electrodes, 'YTickLabel', all_electrodes);
                axis square;
                box off
                set(gca, 'FontSize', 8); % Reduce font size

                hold on;
                %else
                %end

            else
            end

        end
        end


fprintf('Total number of negative correlated pairs:')
size(neg_corrs,1)
%% correlations: baseline only
 

bound1 = 5800; bound2 = 6000;
total_comparisons = 0;
neg_corrs=[];noise_correlation_store = [];null_corrs=[];


    neg_corrs=[];noise_correlation_store = [];all_pairs = [];  % will store ALL tested electrode pairs with p and corr

    for array = 1:2
        % Initialize noise_electrode cell array
        for n = 1:6
            clear I;
            clear coordinate_med;clear coordinate_trans;clear noise_electrode;
            electrodes = find(subject(n).array==array);
            clear noise_electrode;


            % Compute noise correlation matrix and p-values
            num_electrodes = length(electrodes);

            counter = 1;

            for i = 1:num_electrodes
                e = electrodes(i);

                % --- Use all decision 0/1 trials for this subject (same filter as before) ---
which_trials_all = find(subject(n).decision < 2 & ...
    subject(n).electrode(e).trigger(1).good_trials==1);   % baseline trigger only

if isempty(which_trials_all)
    continue;
end

% --- RAW baseline HFA in [bound1:bound2] from trigger(1) ---
BL = subject(n).electrode(e).trigger(1).high_gamma_mat;   % [trials x time]
raw_all = BL(which_trials_all, bound1:bound2);            % [nSelTrials x T]

% If gpuArray, gather
if isa(raw_all,'gpuArray')
    raw_all = gather(raw_all);
end

% Store raw baseline HFA (trials x time)
raw_electrode{counter} = raw_all;

% Coordinates unchanged
coordinate_med(counter)   = subject(n).electrode(e).med_coor;
coordinate_trans(counter) = subject(n).electrode(e).trans_coor;

counter = counter + 1;

            end


            noise_correlation_matrix = zeros(counter-1, counter-1);
            pval_matrix = zeros(counter-1, counter-1);

            if exist('coordinate_med')
                [B,I] = sort(coordinate_med);

                for i = 1:counter-1
                    for j = i:counter-1
                        if i ==j | abs(i - j)<1
                            noise_correlation_matrix(I(i), I(j)) = 1;
                            pval_matrix(I(i), I(j)) = 1;
                        else

                            % Compute Pearson correlation and p-value
                            [corr_value, p_value] = calculate_corr(raw_electrode{I(i)}, raw_electrode{I(j)});

                            % Store in matrices
                            noise_correlation_matrix(i,j) = corr_value;
                            noise_correlation_matrix(j,i) = corr_value; % Symmetric matrix
                            if corr_value < 0
                                neg_corrs = [neg_corrs; n electrodes(I(i)) electrodes(I(j)) B(i) B(j) coordinate_trans(I(i)) coordinate_trans(I(j)) p_value corr_value];
                            else

                            end

                            null_corrs = [null_corrs; n electrodes(I(i)) electrodes(I(j)) B(i) B(j) coordinate_trans(I(i)) coordinate_trans(I(j)) corr_value];
                            all_pairs = [all_pairs; n array electrodes(I(i)) electrodes(I(j)) B(i) B(j) coordinate_trans(I(i)) coordinate_trans(I(j)) p_value corr_value];

                            total_comparisons = total_comparisons + 1;

                            pval_matrix(i,j) = p_value;
                            pval_matrix(j,i) = p_value; % Symmetric matrix
                        end
                    end
                end

                % Convert the matrix into a single column vector
                vector_column1 = noise_correlation_matrix(:);

                % Create a second column with the value of n repeated
                vector_column2 = repmat(n, length(vector_column1), 1);

                % Combine into a two-column matrix
                result_matrix = [vector_column1, vector_column2];

                noise_correlation_store = [noise_correlation_store;result_matrix];

                noise_correlation_matrices{n} = noise_correlation_matrix;
                pval_matrices{n} = log10(pval_matrix);

                % Subplot 1: Noise correlation matrix
                %if n ==4 && array == 2 %plotting example matrix
                subplot(6,3,(n-1)*3+array)
                imagesc(noise_correlation_matrices{n}); % Display correlation matrix
                colorbar;
                colormap gray
                caxis([-0.05 0])
                all_electrodes = electrodes(I);
                num_electrodes = length(all_electrodes);
                set(gca, 'XTick', 1:num_electrodes, 'XTickLabel', all_electrodes);
                set(gca, 'YTick', 1:num_electrodes, 'YTickLabel', all_electrodes);
                axis square;
                box off
                set(gca, 'FontSize', 8); % Reduce font size

                hold on;
                %else
                %end

            else
            end

        end
        end


fprintf('Total number of negative correlated pairs:')
size(neg_corrs,1)
%% =========================
%  Post-hoc reporting block (FIXED)
%  Requires: neg_corrs, total_comparisons
%  Requires for null tests: all_pairs
%  Requires for coeff tests: sorted_coeffs, sorted_subject_ID, sorted_electrode_ID
% ==========================

alpha = 0.05;
alpha_bonf = alpha / max(total_comparisons,1);

fprintf('\n====================\n');
fprintf('Noise-corr summary (bound1=%d, bound2=%d)\n', bound1, bound2);
fprintf('Total correlations computed: %d\n', total_comparisons);

% neg_corrs columns:
% [1 subj, 2 elec1, 3 elec2, 4 med1, 5 med2, 6 trans1, 7 trans2, 8 p, 9 corr]

if isempty(neg_corrs)
    fprintf('No negative correlations found (corr < 0).\n');
    fprintf('====================\n');
    return;
end

p_neg = neg_corrs(:,8);

sig05_idx   = (p_neg < alpha);
sigBonf_idx = (p_neg < alpha_bonf);

fprintf('Negative-correlated pairs (corr < 0): %d\n', size(neg_corrs,1));
fprintf('  Anticorrelated pairs with p < 0.05: %d\n', sum(sig05_idx));
fprintf('  Anticorrelated pairs with Bonferroni p < %.3g: %d\n', alpha_bonf, sum(sigBonf_idx));

% ---------- choose which significance threshold you want to use for A ----------
use_bonf = true;   % <--- toggle this
if use_bonf
    A = neg_corrs(sigBonf_idx,:);
    thresh_str = sprintf('p < %.3g (Bonferroni)', alpha_bonf);
else
    A = neg_corrs(sig05_idx,:);
    thresh_str = 'p < 0.05';
end

if isempty(A)
    fprintf('\nNo anticorrelated pairs at %s.\n', thresh_str);
    fprintf('====================\n');
    return;
end

%% ---------- Spatial summary for anticorrelated pairs ----------
med1 = A(:,4); med2 = A(:,5);
med_moreMedial  = min(med1, med2);
med_moreLateral = max(med1, med2);

fprintf('\nSpatial summary for anticorrelated pairs (%s):\n', thresh_str);
fprintf('  Mean coordinate of MORE medial electrode:  %.3f\n', mean(med_moreMedial,'omitnan'));
fprintf('  Mean coordinate of MORE lateral electrode: %.3f\n', mean(med_moreLateral,'omitnan'));

%% ==========================================================
%  Null comparison: ML separation of anticorr vs ALL pairs
% ==========================================================

if ~exist('all_pairs','var') || isempty(all_pairs)
    error('all_pairs not found. Cannot construct null distribution.');
end

% Anticorr ML separation
ML_sep_anticorr = max(A(:,4:5),[],2) - min(A(:,4:5),[],2);
ML_sep_anticorr = ML_sep_anticorr(~isnan(ML_sep_anticorr));
mean_anticorr = mean(ML_sep_anticorr,'omitnan');

% All-pairs ML separation
ML_sep_all = max(all_pairs(:,5:6),[],2) - min(all_pairs(:,5:6),[],2);
ML_sep_all = ML_sep_all(~isnan(ML_sep_all));
mean_null = mean(ML_sep_all,'omitnan');

% Global permutation null
Nperm = 10000;
nA = numel(ML_sep_anticorr);
perm_means = nan(Nperm,1);
for p = 1:Nperm
    idx = randsample(numel(ML_sep_all), nA, true);
    perm_means(p) = mean(ML_sep_all(idx),'omitnan');
end
p_perm_global = mean(perm_means >= mean_anticorr);

fprintf('\nSpatial segregation test (medial-lateral):\n');
fprintf('  Mean ML separation (anticorr pairs): %.3f\n', mean_anticorr);
fprintf('  Mean ML separation (all pairs):      %.3f\n', mean_null);
fprintf('  Permutation p-value (anticorr > null): %.4f\n', p_perm_global);

%% ---------- Subject-matched permutation null ----------
subjects = unique(A(:,1));
perm_means = nan(Nperm,1);

for p = 1:Nperm
    temp = [];
    for s = subjects'
        idxA  = (A(:,1) == s);
        n_s   = sum(idxA);
        if n_s == 0, continue; end

        idxAll = (all_pairs(:,1) == s);
        if ~any(idxAll), continue; end

        ML_s = max(all_pairs(idxAll,5:6),[],2) - min(all_pairs(idxAll,5:6),[],2);
        ML_s = ML_s(~isnan(ML_s));
        if isempty(ML_s), continue; end

        temp = [temp; ML_s(randsample(numel(ML_s), n_s, true))];
    end
    perm_means(p) = mean(temp,'omitnan');
end

p_perm_subject = mean(perm_means >= mean_anticorr);
fprintf('  Subject-matched permutation p-value: %.4f\n', p_perm_subject);

%% ==========================================================
%  Reverse conditioning (UNIQUE PAIRS only)
%  P(more medial | coeff > 0) and P(more lateral | coeff < 0)
% ==========================================================

if exist('sorted_coeffs','var') && exist('sorted_subject_ID','var') && exist('sorted_electrode_ID','var')

    subjA = A(:,1);
    e1    = A(:,2);
    e2    = A(:,3);
    med1A = A(:,4);
    med2A = A(:,5);

    % Unique by (subject, unordered electrode pair)
    emin = min(e1,e2);
    emax = max(e1,e2);
    [~, ia] = unique([subjA emin emax], 'rows', 'stable');

    A_u    = A(ia,:);
    subj_u = A_u(:,1);
    e1_u   = A_u(:,2);
    e2_u   = A_u(:,3);
    med1_u = A_u(:,4);
    med2_u = A_u(:,5);

    % medial vs lateral per unique pair
    medial_elec  = e1_u;
    lateral_elec = e2_u;
    swap = med2_u < med1_u;          % e2 is more medial
    medial_elec(swap)  = e2_u(swap);
    lateral_elec(swap) = e1_u(swap);

    % coeff lookup
    coeff_med = nan(size(medial_elec));
    coeff_lat = nan(size(lateral_elec));

    for k = 1:numel(medial_elec)
        s = subj_u(k);

        idxm = find(sorted_subject_ID == s & sorted_electrode_ID == medial_elec(k), 1, 'first');
        if ~isempty(idxm), coeff_med(k) = sorted_coeffs(idxm); end

        idxl = find(sorted_subject_ID == s & sorted_electrode_ID == lateral_elec(k), 1, 'first');
        if ~isempty(idxl), coeff_lat(k) = sorted_coeffs(idxl); end
    end

    % electrode-role table (2*#pairs)
    role_is_medial  = [true(numel(coeff_med),1);  false(numel(coeff_lat),1)];
    role_is_lateral = [false(numel(coeff_med),1); true(numel(coeff_lat),1)];

    coeff_all = [coeff_med(:); coeff_lat(:)];
    valid_all = ~isnan(coeff_all);

    idx_pos = valid_all & (coeff_all > 0);
    idx_neg = valid_all & (coeff_all < 0);

    prop_med_given_pos = 100 * sum(role_is_medial(idx_pos)) / max(sum(idx_pos),1);
    prop_lat_given_neg = 100 * sum(role_is_lateral(idx_neg)) / max(sum(idx_neg),1);

    fprintf('\nReverse conditioning within anticorr pairs (%s), UNIQUE PAIRS:\n', thresh_str);
    fprintf('  Unique anticorr pairs used: %d (from %d rows in A)\n', numel(ia), size(A,1));
    fprintf('  Of electrodes with sorted_coeff > 0, %% that are MORE MEDIAL: %.1f%% (%d/%d)\n', ...
        prop_med_given_pos, sum(role_is_medial(idx_pos)), sum(idx_pos));
    fprintf('  Of electrodes with sorted_coeff < 0, %% that are MORE LATERAL: %.1f%% (%d/%d)\n', ...
        prop_lat_given_neg, sum(role_is_lateral(idx_neg)), sum(idx_neg));

else
    fprintf('\nNOTE: sorted_coeffs / sorted_subject_ID / sorted_electrode_ID not found.\n');
    fprintf('      Skipping coefficient sign analysis.\n');
end

fprintf('====================\n');


%% heatmap of spatial distribution of the more medially located electrode in each negatively correlated pair

% Parameters
window_size = 10;  % Size of the sliding window in mm
overlap = 5;       % Overlap between windows in mm

% Extract X and Y coordinates
X_coords = neg_corrs(:, 4);
Y_coords = neg_corrs(:, 6);

% minimum of ranges covered by electrodes from all subjects +/- overlap
x_min = -23;
x_max = 43;
y_min = -23;
y_max = 30;

% Generate sliding window centers with overlap
x_centers = x_min:overlap:x_max;
y_centers = y_min:overlap:y_max;

% Initialize matrix to store electrode counts in each bin
electrode_counts = zeros(length(y_centers), length(x_centers));

% Extract X and Y coordinates
X_coords = neg_corrs(neg_corrs(:,8)<sigBonf_idx, 4);
Y_coords = neg_corrs(neg_corrs(:,8)<sigBonf_idx, 6);

% Loop over sliding window centers to calculate the electrode count in each bin
for x_idx = 1:length(x_centers)
    for y_idx = 1:length(y_centers)
        % Define the window range
        x_start = x_centers(x_idx) - window_size/2;
        x_end = x_centers(x_idx) + window_size/2;
        y_start = y_centers(y_idx) - window_size/2;
        y_end = y_centers(y_idx) + window_size/2;

        % Count electrodes within the current window
        in_window = X_coords >= x_start & X_coords < x_end & ...
                    Y_coords >= y_start & Y_coords < y_end;
        electrode_counts(y_idx, x_idx) = sum(in_window);
    end
end

% Set bins with no electrodes to NaN for white background
plot_data = electrode_counts;
plot_data(electrode_counts == 0) = NaN;  % Set empty bins to NaN

% Plot the heatmap
figure;
h = imagesc(x_centers, y_centers, plot_data);  % Plot the data
colormap([0 0 0 ]);  % Set hot pink as the colormap
%colorbar off;  % Turn off the colorbar since color doesn't vary

% Set transparency based on the number of electrodes in each bin
alpha_data = electrode_counts ./ max(max(electrode_counts)); % Normalize alpha
set(h, 'AlphaData', alpha_data);

% Set axis limits
xlim([-25, 45]);
ylim([-25, 30]);

% Ensure y-axis is not flipped
set(gca, 'YDir', 'normal');

% Add labels and title
title('Electrode Density Heatmap', 'FontSize', 14);
xlabel('X Coordinate (mm)', 'Fontsize', 14);
ylabel('Y Coordinate (mm)', 'Fontsize', 14);
set(gca, 'TickDir', 'out');
box off;
hold on;
%% heatmap of spatial distribution of the more laterally located electrode in each negatively correlated pair
figure(3)
% Parameters
window_size = 10;  % Size of the sliding window in mm
overlap = 5;       % Overlap between windows in mm

% Extract X and Y coordinates
X_coords = neg_corrs(:, 4);
Y_coords = neg_corrs(:, 6);

% Determine the min and max coordinates for the grid
x_min = floor(min(X_coords)) - window_size;
x_max = ceil(max(X_coords)) + window_size;
y_min = floor(min(Y_coords)) - window_size;
y_max = ceil(max(Y_coords)) + window_size;

% Generate sliding window centers with overlap
x_centers = x_min:overlap:x_max;
y_centers = y_min:overlap:y_max;

% Initialize matrix to store electrode counts in each bin
electrode_counts = zeros(length(y_centers), length(x_centers));

% Extract X and Y coordinates
X_coords = neg_corrs(neg_corrs(:,8)<sigBonf_idx, 5);
Y_coords = neg_corrs(neg_corrs(:,8)<sigBonf_idx, 7);

% Loop over sliding window centers to calculate the electrode count in each bin
for x_idx = 1:length(x_centers)
    for y_idx = 1:length(y_centers)
        % Define the window range
        x_start = x_centers(x_idx) - window_size/2;
        x_end = x_centers(x_idx) + window_size/2;
        y_start = y_centers(y_idx) - window_size/2;
        y_end = y_centers(y_idx) + window_size/2;

        % Count electrodes within the current window
        in_window = X_coords >= x_start & X_coords < x_end & ...
                    Y_coords >= y_start & Y_coords < y_end;
        electrode_counts(y_idx, x_idx) = sum(in_window);
    end
end

% Set bins with no electrodes to NaN for white background
plot_data = electrode_counts;
plot_data(electrode_counts == 0) = NaN;  % Set empty bins to NaN

% Plot the heatmap
figure;
h = imagesc(x_centers, y_centers, plot_data);  % Plot the data
colormap([0 0 0 ]);  % Set hot pink as the colormap
%colorbar off;  % Turn off the colorbar since color doesn't vary

% Set transparency based on the number of electrodes in each bin
alpha_data = electrode_counts ./ max(max(electrode_counts)); % Normalize alpha
set(h, 'AlphaData', alpha_data);

% Set axis limits
xlim([-25, 45]);
ylim([-25, 30]);

% Ensure y-axis is not flipped
set(gca, 'YDir', 'normal');

% Add labels and title
title('Electrode Density Heatmap', 'FontSize', 14);
xlabel('X Coordinate (mm)', 'Fontsize', 14);
ylabel('Y Coordinate (mm)', 'Fontsize', 14);
set(gca, 'TickDir', 'out');
box off;
hold on;


%%
%% ==========================================================
%  Lag analysis with padded time window
% ==========================================================

L = 400;                      % max lag in samples (ms)
lags = -L:10:L;

core_idx = bound1:bound2;     % e.g., 5700:6000
ext_idx  = (bound1-L):(bound2+L);

nCore = numel(core_idx);

corr_by_pair = nan(size(A,1), numel(lags));

for k = 1:size(A,1)

    s  = A(k,1);
    e1 = A(k,2);
    e2 = A(k,3);

    % --- matched trials ---
    which1 = find(subject(s).decision < 2 & ...
                  subject(s).electrode(e1).trigger(3).good_trials==1 & ...
                  subject(s).electrode(e1).trigger(1).good_trials==1);

    which2 = find(subject(s).decision < 2 & ...
                  subject(s).electrode(e2).trigger(3).good_trials==1 & ...
                  subject(s).electrode(e2).trigger(1).good_trials==1);

    which = intersect(which1, which2);
    if isempty(which), continue; end

    % --- value + decision labels ---
    ttt = subject(s).trial_type_trial;
    vtt = subject(s).value_trial_type;
    val_per_trial = vtt(ttt);

    dec = subject(s).decision(which);
    val = val_per_trial(which);

    % --- pull matrices ---
    HG1 = subject(s).electrode(e1).trigger(3).high_gamma_mat;
    BL1 = subject(s).electrode(e1).trigger(1).high_gamma_mat;
    HG2 = subject(s).electrode(e2).trigger(3).high_gamma_mat;
    BL2 = subject(s).electrode(e2).trigger(1).high_gamma_mat;

    % Make sure ext_idx is within available time range
    T1 = size(HG1,2);
    T2 = size(HG2,2);
    ext_idx1 = ext_idx(ext_idx>=1 & ext_idx<=T1);
    ext_idx2 = ext_idx(ext_idx>=1 & ext_idx<=T2);

    % If your data don't cover the full padded range, skip this pair
    if numel(ext_idx1) < numel(ext_idx) || numel(ext_idx2) < numel(ext_idx)
        continue;
    end

    % --- baseline scalar per trial (same as before) ---
    b1 = mean(BL1(which, 5500:6000), 2);
    b2 = mean(BL2(which, 5500:6000), 2);

    % --- baseline-normalized HFA on EXTENDED window ---
    sig1_ext = HG1(which, ext_idx);
    sig2_ext = HG2(which, ext_idx);

    if isa(sig1_ext,'gpuArray'), sig1_ext = gather(sig1_ext); end
    if isa(sig2_ext,'gpuArray'), sig2_ext = gather(sig2_ext); end
    if isa(b1,'gpuArray'), b1 = gather(b1); end
    if isa(b2,'gpuArray'), b2 = gather(b2); end

    norm1_ext = (sig1_ext - b1) ./ b1;   % [trials x extT]
    norm2_ext = (sig2_ext - b2) ./ b2;

    % --- residualize by (decision,value) condition on EXT window ---
    cond_mat = [dec(:), val(:)];
    [~,~,cond_id] = unique(cond_mat,'rows');
    nConds = max(cond_id);

    noise1_ext = zeros(size(norm1_ext));
    noise2_ext = zeros(size(norm2_ext));

    for c = 1:nConds
        idxc = (cond_id == c);
        noise1_ext(idxc,:) = norm1_ext(idxc,:) - mean(norm1_ext(idxc,:),1);
        noise2_ext(idxc,:) = norm2_ext(idxc,:) - mean(norm2_ext(idxc,:),1);
    end

    % Pre-extract electrode 2 CORE window from extended indices
    % Map core indices into ext_idx coordinates
    core_in_ext = (core_idx - (bound1-L) + 1);   % position inside ext window (1..extLen)
    E2_core = noise2_ext(:, core_in_ext);        % [trials x nCore]

    % --- lag sweep: shift ELECTRODE 1 window relative to electrode 2 core ---
    for li = 1:numel(lags)
        lag = lags(li);

        % Electrode 1 window shifted by lag, still length nCore
        idx_shift = core_in_ext + lag;

        % Because we padded by ±L, idx_shift stays in-bounds for all lag in [-L,L]
        E1_shift = noise1_ext(:, idx_shift);

        % Correlate vectorized trial×time samples
        x = E1_shift(:);
        y = E2_core(:);
        good = ~isnan(x) & ~isnan(y);

        if nnz(good) > 10
            corr_by_pair(k,li) = corr(x(good), y(good), 'type','Pearson');
        end
    end
end

% summarize and plot
mean_r = mean(corr_by_pair, 1, 'omitnan');
sem_r  = std(corr_by_pair, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(corr_by_pair),1));
%% ==========================================================
%  Summarize, pretty plot (patch), and compute offset statistics
% ==========================================================

% summarize across pairs at each lag
mean_r = mean(corr_by_pair, 1, 'omitnan');

n_eff = sum(~isnan(corr_by_pair),1);                   % how many pairs contribute at each lag
sem_r = std(corr_by_pair, 0, 1, 'omitnan') ./ sqrt(max(n_eff,1));

% ---- Pretty plot: patch = mean ± SEM ----
figure; hold on;

x = lags(:);
y = mean_r(:);
e = sem_r(:);

% patch polygon
xp = [x; flipud(x)];
yp = [y - e; flipud(y + e)];
hpatch = patch(xp, yp, 0.8*[1 1 1], 'EdgeColor','none', 'FaceAlpha',0.4); %#ok<NASGU>

% mean line
plot(x, y, 'k-', 'LineWidth', 2);

% zero line
yline(0,'k-');

xlabel('Temporal offset applied to electrode 1 (ms)');
ylabel('Correlation (r)');
title(sprintf('Lag dependence of anticorrelation (pairs with any data, N=%d)', ...
    sum(any(~isnan(corr_by_pair),2))));

set(gca,'TickDir','out');
box off;

%% ----------------------------------------------------------
%  Stats: compare |lag|<100ms vs outer ranges
%  (paired across electrode pairs)
% ----------------------------------------------------------

% Define lag bins
idx_center = (lags >= -100 & lags <= 100);
idx_left   = (lags >= -400 & lags <  -200);
idx_right  = (lags >   200 & lags <= 400);

% Pairwise averages within each lag-range (one value per pair)
rho_center = mean(corr_by_pair(:, idx_center), 2, 'omitnan');
rho_left   = mean(corr_by_pair(:, idx_left),   2, 'omitnan');
rho_right  = mean(corr_by_pair(:, idx_right),  2, 'omitnan');

% Only keep pairs with both values present for each comparison
vL = ~isnan(rho_center) & ~isnan(rho_left);
vR = ~isnan(rho_center) & ~isnan(rho_right);

% Paired t-tests (center vs left, center vs right)
% Your text implies rho_center is "smaller" (more negative) than outer ranges.
if any(vL)
    [~,pL,~,stL] = ttest(rho_center(vL), rho_left(vL));  % tests mean(center-left)=0
    fprintf('\nCenter vs Left:\n');
    fprintf('  N pairs: %d\n', sum(vL));
    fprintf('  mean rho center: %.4f\n', mean(rho_center(vL),'omitnan'));
    fprintf('  mean rho left:   %.4f\n', mean(rho_left(vL),'omitnan'));
    fprintf('  t(%d)=%.2f, p=%.4g\n', stL.df, stL.tstat, pL);
else
    fprintf('\nCenter vs Left: not enough data.\n');
end

if any(vR)
    [~,pR,~,stR] = ttest(rho_center(vR), rho_right(vR));
    fprintf('\nCenter vs Right:\n');
    fprintf('  N pairs: %d\n', sum(vR));
    fprintf('  mean rho center: %.4f\n', mean(rho_center(vR),'omitnan'));
    fprintf('  mean rho right:  %.4f\n', mean(rho_right(vR),'omitnan'));
    fprintf('  t(%d)=%.2f, p=%.4g\n', stR.df, stR.tstat, pR);
else
    fprintf('\nCenter vs Right: not enough data.\n');
end

