%% Build Z-scored heatmap of HFA with Optimized GPU Usage
folders = {'EC304' 'EC288' 'PR05' 'PR06' 'BJH058' 'DP01'};
% Preallocate GPU memory for large variables
num_electrodes=[];
for i=1:length(folders)
    num_electrodes = [num_electrodes length(subject(i).electrode)];
end

bin_width = 20; % User-specified bin width in milliseconds
num_bins = 3000 / bin_width ;  % Approximate number of bins
all_z_scores = zeros(sum(num_electrodes), num_bins);  % Keep this on CPU
brain_region = [];

all_subjects = {};
electrode_ID = [];
subject_ID = [];
olf_coor=[];lat_coor = [];mid_coor=[];med_coor = [];trans_coor = [];
side = [];

t = 3;

% Loop through each subject
for n = 1:length(folders)
    num_electrodes_subject = num_electrodes(n);

    % Loop through each electrode
    for e = 1:num_electrodes_subject

        % Move entire high_gamma_mat to GPU once (avoiding CPU-GPU bottleneck)
        signal_gpu = gpuArray(subject(n).electrode(e).trigger(t).high_gamma_mat(:, 4500:7499));

        % exclude any skipped trials
        trialIDs1 = find(subject(n).decision<2 & ~isnan(signal_gpu(:,1)));

        % GPU-based binning
        binned_data1_gpu = bin_data_columns(bin_width, trialIDs1, signal_gpu);

        % Process baseline signal (GPU)
        signal_baseline_gpu = gpuArray(subject(n).electrode(e).trigger(1).high_gamma_mat(:, 5700:5999));
        basedistr_gpu = bin_data_columns(bin_width, trialIDs1, signal_baseline_gpu);

        % After you gather to CPU:
        binned_data_cpu = gather(binned_data1_gpu);   % [nTrials x nSignalBins]
        basedistr_cpu   = gather(basedistr_gpu);      % [nTrials x nBaseBins]

        % One baseline number per trial (mean over the baseline window)
        baseline_per_trial = mean(basedistr_cpu, 2);  % [nTrials x 1]

        ntrials = size(binned_data_cpu, 1);
        z_tc = zeros(1, size(binned_data_cpu, 2));

        for j = 1:size(binned_data_cpu, 2)
            signal_per_trial = binned_data_cpu(:, j);                  % [nTrials x 1]
            d = signal_per_trial - baseline_per_trial;                 % paired diffs
            % z of the mean difference (normal approx to t):
            sd_d = std(d, 0);                                         % sample SD of diffs
            z_tc(j) = mean(d) / (sd_d / sqrt(ntrials ));
        end


        % Store results in CPU array (no need for GPU at this step)
        if n ==1
        all_z_scores(e, :) = z_tc;
        else
            all_z_scores(e + sum(num_electrodes(1:n-1)), :) = z_tc;
        end

        %side = [side; subject(n).electrode(e).brain_region(1)=='R'];
        brain_region = [brain_region; {subject(n).electrode(e).anatomical_info}];
        olf_coor = [olf_coor; subject(n).electrode(e).olf_coor];
        lat_coor = [lat_coor; subject(n).electrode(e).lat_coor];
        med_coor = [med_coor; subject(n).electrode(e).med_coor];
        trans_coor = [trans_coor; subject(n).electrode(e).trans_coor];
        mid_coor = [mid_coor; subject(n).electrode(e).mid_coor];

        % Keep track of metadata
        all_subjects{end+1} = folders{n};
        subject_ID(end+1) = n;
        electrode_ID(end+1) = e;

        % Clear GPU memory to prevent memory overload
        clear signal_gpu signal_baseline_gpu binned_data1_gpu basedistr_gpu;
    end
end

%% GLM for individual electrode encoding of decision
if exist('coeffs')
    clear coeffs; clear pvals;
end


for j=1:length(subject_ID)
    e = electrode_ID(j);n = subject_ID(j);
    filter = find(subject(n).decision<2  & subject(n).electrode(e).trigger(3).good_trials==1 & subject(n).electrode(e).trigger(1).good_trials==1);
    signal = [];response = [];baseline=[];
    for a = 1:length(filter)
        signal = [signal; sum(subject(n).electrode(e).trigger(3).high_gamma_mat(filter(a),5700:6000),2)/(301)];
        baseline = [baseline;sum(subject(n).electrode(e).trigger(1).high_gamma_mat(filter(a),5500:6000),2)/301];

        response = [response; subject(n).decision(filter(a))...
        subject(n).value_trial(filter(a))]; %<--to examine decision
        %encoding independent of subjective value
    end

    % Create the GLM
    glm_model = fitglm(response,(signal-baseline)./baseline,'linear');

    num_coeffs = size(glm_model.Coefficients)-1;
    for b = 1:num_coeffs
        coeffs{b}(j) = glm_model.Coefficients.Estimate(b+1);
        pvals{b}(j) = glm_model.Coefficients.pValue(b+1);
    end
end

%% generate flags for excluding any electrodes

% Exclude certain indices
% Convert cell array to logical mask (1 if 'posterior OFC', 0 otherwise)
mask = olf_coor>0;
exclude_indices = [find(~mask)];

% Create new variables to store the flagged data
flagged_z_scores = all_z_scores;
flagged_subjects = all_subjects;
flagged_subject_ID = subject_ID;
flagged_electrode_ID = electrode_ID;
flagged_brain_region = brain_region;
flagged_olf_coor = olf_coor;
flagged_lat_coor = lat_coor;
flagged_med_coor = med_coor;
flagged_trans_coor = trans_coor;
flagged_mid_coor = mid_coor;
flagged_coeffs = coeffs{1};

% Remove the excluded indices from flagged_roc_tc and flagged_coordinates
flagged_z_scores(exclude_indices, :) = [];
flagged_brain_region(exclude_indices) = [];
flagged_subjects(exclude_indices) = [];
flagged_subject_ID(exclude_indices) = [];
flagged_electrode_ID(exclude_indices) = [];
flagged_olf_coor(exclude_indices) = [];
flagged_mid_coor(exclude_indices) = [];
flagged_lat_coor(exclude_indices) = [];
flagged_med_coor(exclude_indices) = [];
flagged_trans_coor(exclude_indices) = [];
flagged_coeffs(exclude_indices) = [];



%% sort by anatomical order
[~, plotorder] = sort(flagged_med_coor);
sorted_z_scores = flagged_z_scores(plotorder, :);
sorted_subjects = flagged_subjects(plotorder);
sorted_subject_ID = flagged_subject_ID(plotorder);
sorted_electrode_ID = flagged_electrode_ID(plotorder);
sorted_med_coor = flagged_med_coor(plotorder);
sorted_lat_coor = flagged_lat_coor(plotorder);
sorted_trans_coor = flagged_trans_coor(plotorder);
sorted_mid_coor = flagged_mid_coor(plotorder);
sorted_coeffs = flagged_coeffs(plotorder);

%% plot heatmap of z scores sorted anatomically

% --- Subject grayscale colors ---
for i = 1:n
    subject_colors(i,:) = [1-i/n, 1-i/n, 1-i/n];
end

figure(1); clf
standard_font_size = 32; numbers_font_size = 20;

mu = mean(sorted_z_scores(:));
sd = std(sorted_z_scores(:));
upper = mu + 2*sd;
lower = mu - 2*sd;

% ================= HEATMAP (keeps its own colormap) =================
ax_heat = axes('Position', [0.05, 0.1, 0.55, 0.8]);
imagesc(ax_heat, sorted_z_scores);
hold(ax_heat, 'on');
clim(ax_heat, [lower upper]);
title(ax_heat, 'Z-score values', 'FontSize', standard_font_size);
xlabel(ax_heat, 'Time (bins)', 'FontSize', standard_font_size);
ylabel(ax_heat, 'Electrodes',   'FontSize', standard_font_size);
set(ax_heat, 'FontSize', numbers_font_size, 'TickDir', 'out');
colormap(ax_heat, bone);             % <-- heatmap stays bone (not pink/yellow)
colormap(ax_heat, bone);  % heatmap stays bone
cb = colorbar(ax_heat, 'EastOutside');   % or 'SouthOutside'
cb.Label.String = 'z-score';

% ================= SUBJECT COLOR BLOCK (truecolor) ===================
subject_color_block = zeros(length(sorted_subjects), 1, 3);
for i = 1:length(folders)
    idx = find(strcmp(sorted_subjects, folders{i}));  % rows for this subject
    % Fix the size-mismatch by expanding the 1x1x3 color across all rows in idx
    color13 = reshape(subject_colors(i,:), [1 1 3]);         % 1x1x3
    subject_color_block(idx,1,:) = repmat(color13, [numel(idx) 1 1]);  % Nx1x3
end

ax_subj = axes('Position', [0.62, 0.1, 0.03, 0.8]);          % narrower strip
image(ax_subj, subject_color_block);                         % truecolor image
set(ax_subj, 'YDir', get(ax_heat,'YDir'));                   % align rows
title(ax_subj, 'Subjects', 'FontSize', standard_font_size);
set(ax_subj, 'XTick', [], 'YTick', [], 'TickDir','out', 'FontSize', numbers_font_size);
box(ax_subj, 'off');

% ================= MEDIAL-COORDINATE COLOR STRIP (pink->yellow) =====
% Build pink->yellow colormap used ONLY for the med-coord strip
ncolors = 256;
pink   = [1, 0.2, 0.6];
yellow = [1, 1, 0];
cmap_py = [linspace(pink(1),yellow(1),ncolors)', ...
           linspace(pink(2),yellow(2),ncolors)', ...
           linspace(pink(3),yellow(3),ncolors)'];

vals = sorted_med_coor(:);
vmin = min(vals);
vmax = max(vals);
idxc = max(1, min(ncolors, round( rescale(vals, 1, ncolors) )));
rgb  = cmap_py(idxc, :);                               % num_rows x 3
med_color_block = reshape(rgb, [numel(vals) 1 3]);     % num_rows x 1 x 3

% Vertical strip immediately to the right of subject block
ax_med = axes('Position', [0.66, 0.1, 0.03, 0.8]);
image(ax_med, med_color_block);                        % truecolor, but we still set caxis/colormap for the colorbar scale
set(ax_med, 'YDir', get(ax_heat,'YDir'));              % align rows
set(ax_med, 'XTick', [], 'TickDir', 'out', 'FontSize', numbers_font_size);
ylabel(ax_med, 'Medial (mm)', 'FontSize', standard_font_size);
box(ax_med, 'off');

% Apply the pink->yellow colormap & range to THIS axis for colorbar scaling
colormap(ax_med, cmap_py);
caxis(ax_med, [vmin vmax]);

% ---- Label only 5 mm steps on this strip ----
num_rows = size(sorted_z_scores, 1);
min_mm = floor(vmin/5)*5;
max_mm = ceil(vmax/5)*5;
targets = min_mm:5:max_mm;
[~, nearest_rows] = arrayfun(@(t) min(abs(vals - t)), targets);
[rows_unique, ia] = unique(nearest_rows, 'stable');
targets_unique = targets(ia);
set(ax_med, 'YTick', rows_unique, 'YTickLabel', string(targets_unique), 'YAxisLocation', 'left');
set(ax_med, 'TickLength', [0.01 0.01]);

% Vertical colorbar next to the med-coord strip (RIGHT side)
cb_med = colorbar(ax_med, 'EastOutside');
cb_med.Label.String = 'Medial (mm)';
% Optional: tighten the colorbar placement
cb_med.Position = [0.70, 0.1, 0.02, 0.8];   % [x y w h], adjust as needed

%% plot coefficients of decision encoding with sorted coefficients
x = 1:numel(sorted_coeffs);
y = smooth(sorted_coeffs, 20);  % moving average window = 10

figure; hold on;

% Split into positive and negative parts relative to zero
y_pos = max(y, 0);    % portions above 0
y_neg = min(y, 0);    % portions below 0

% (Optional) avoid NaN issues in area filling
y_pos(isnan(y_pos)) = 0;
y_neg(isnan(y_neg)) = 0;

% Fill areas
h1 = area(x, y_pos, 'BaseValue', 0);
set(h1, 'FaceColor', 'b', 'EdgeColor', 'none', 'FaceAlpha', 0.5);

h2 = area(x, y_neg, 'BaseValue', 0);
set(h2, 'FaceColor', 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.5);

% Plot the smoothed line and the zero line
plot(x, y, 'k', 'LineWidth', 1.5);
yline(0, 'k-');

ylim([-0.06 0.06])


xlabel('Index');
ylabel('Smoothed coeff');
title('Smoothed coefficients with positive/negative fill');
box off; hold off;

%% look at group PSTH
clear trialset
figure(1)
clust_indices = find(abs(sorted_med_coor)<5); %MOS adjacent
%clust_indices = find(sorted_med_coor>(round(max(sorted_med_coor))-15)); %lat-most 1.5cm

t = 3;num_conditions = 2;
signal_bounds = [4000 8000];
baseline_bounds = [5500 6000];
window = signal_bounds(2) - signal_bounds(1) + 1;
smoothing_kernel = 1;
x2 = round(signal_bounds(1):smoothing_kernel:signal_bounds(2));

red_line = [230, 41, 41]/255;
red_fill = [247, 189, 189]/255;

blue_line = [18, 6, 255]/255;
blue_fill = [180, 176, 255]/255;

for j = 1:num_conditions

    if j == 1
        line_color = red_line;
        fill_color = red_fill;
    else

        line_color = blue_line;
        fill_color = blue_fill;
    end

    % fill_color = [j /num_conditions 1 - j /num_conditions 1];
    % line_color = [j /num_conditions 1 - j /num_conditions 1]*0.6;

    for_psth = zeros(length(clust_indices), window);

    for i = 1:length(clust_indices)

         n = sorted_subject_ID(clust_indices(i));
         e = sorted_electrode_ID(clust_indices(i)); 

        % **IMPORTANT: good_trials filters out noisy trials, or those where
        % epileptogenic activity was noted in pre-processing

        good =  subject(n).electrode(e).trigger(t).good_trials==1 & subject(n).electrode(e).trigger(1).good_trials==1;

        trialset = [{subject(n).decision==0  & good}, ...
            {subject(n).decision==1 & good}];

        trialIDs=find(trialset{j});

            % 'temp' gathers what the electrode does for all trials
            temp = [subject(n).electrode(e).trigger(t).high_gamma_mat(trialIDs, signal_bounds(1):signal_bounds(2))];

            % baseline
            temp_baseline = [subject(n).electrode(e).trigger(1).high_gamma_mat(trialIDs, baseline_bounds(1):baseline_bounds(2))];
            baseline = (sum(sum(temp_baseline))/(size(temp_baseline,1)* size(temp_baseline,2)));

            for_psth(i, :) = (sum(temp,1)/size(temp,1)-baseline)/baseline;

    end

    for_psth = for_psth(~isnan(for_psth(:,1)),:);
    y = (sum(for_psth) / size(for_psth, 1));
    dy = zeros(1,length(y));
    for c = 1:size(for_psth,2)
        dy(c) = std(for_psth(:, c))/sqrt(size(for_psth,1));
    end

    for_anova{j} = for_psth(~isnan(for_psth(:,1)),:);

        fill([(signal_bounds(1):signal_bounds(2))-6000, fliplr((signal_bounds(1):signal_bounds(2))-6000)], [(y - dy)  fliplr((y + dy))],fill_color, 'linestyle', 'none','FaceAlpha', 0.6);
        hold on
        plot((signal_bounds(1):signal_bounds(2))'-6000, y, 'Linewidth',2,'Color',line_color)
        hold on


end

% Adjust plot appearance
set(gca, 'TickDir', 'out');
% Get the current x-axis limits
x_limits = xlim;


% Get the current y-axis limits
y_limits = ylim;


box off;
set(gcf, 'Color', 'w');
xlabel('time from decision (milliseconds)', 'FontSize', 20);
ylabel('\Delta HFA', 'Interpreter', 'tex', 'FontSize', 20);

xline(0,'k--','Linewidth',2)
xlim([-1500 500])

[h p c stats]=ttest([sum(for_anova{1}(:,1700:2000),2)'],[sum(for_anova{2}(:,1700:2000),2)'])


%% scatterplot medial_coordinate versus beta_decision
x = sorted_med_coor;
y = sorted_coeffs';

figure; hold on;
plot(x, y, 'k.', 'MarkerSize', 10);   % points

% Best-fit (least-squares) line
p = polyfit(x, y, 1);
xx = linspace(min(x), max(x), 200);
yy = polyval(p, xx);
plot(xx, yy, 'r-', 'LineWidth', 2);

% --- Correlations ---
% Pearson
[Rp, Pp] = corrcoef(x, y);
r_pearson = Rp(1,2);
p_pearson = Pp(1,2);

% Spearman
[r_spear, p_spear] = corr(x, y, 'Type', 'Spearman');

% --- Style ---
box off;
set(gca, 'TickDir', 'out');
set(gcf, 'Color', 'w');

% --- Annotate (top-left) ---
xl = xlim; yl = ylim;
txt = sprintf(['Pearson r = %.2f, p = %.3g\n' ...
               'Spearman \\rho = %.2f, p = %.3g'], ...
               r_pearson, p_pearson, r_spear, p_spear);

text(xl(1), yl(2), txt, 'VerticalAlignment','top', ...
     'HorizontalAlignment','left', 'FontSize', 12);

xlabel('Medial/Lateral coordinate');
ylabel('\beta (decision ~ activity)');
title('\beta vs. medial/lateral coordinate');

%% scatterplot medial_coordinate versus z-score just prior to a decision
% note: need to re-run the z-scoring analysis with 300ms chunks instead of
% 20ms to reproduce the exact analysis in the manuscript figure 2e
x = sorted_med_coor;
y = sorted_z_scores(:,5);

figure; hold on;

% Scatter points
plot(x, y, 'k.', 'MarkerSize', 10);

% Best-fit (least-squares) line
p = polyfit(x, y, 1);
xx = linspace(min(x), max(x), 200);
yy = polyval(p, xx);
plot(xx, yy, 'r-', 'LineWidth', 2);

% --- Correlations ---
% Pearson
[Rp, Pp] = corrcoef(x, y);
r_pearson = Rp(1,2)
p_pearson = Pp(1,2)

% Spearman
[r_spear, p_spear] = corr(x, y, 'Type', 'Spearman')

% --- Style ---
box off;
set(gca, 'TickDir', 'out');
set(gcf, 'Color', 'w');

% --- Annotate results ---
xl = xlim; yl = ylim;

txt = sprintf(['Pearson r = %.2f, p = %.3g\n' ...
               'Spearman ρ = %.2f, p = %.3g'], ...
               r_pearson, p_pearson, r_spear, p_spear);

text(xl(1), yl(2), txt, 'VerticalAlignment','top', ...
     'HorizontalAlignment','left', 'FontSize', 12);

xlabel('Medial/Lateral coordinate');
ylabel('z-score prior to decision');
title('z-score prior to decision vs. medial/lateral coordinate');


%% look at individual channel PSTH
n = 4; % Subject index
e = 8; % 

t = 3; % Trigger index
good =  subject(n).electrode(e).trigger(t).good_trials==1 & subject(n).electrode(e).trigger(1).good_trials==1;

% Define signal and baseline bounds
signal_bounds = [3000 9000];
baseline_bounds = [5500 6000];


% Define trial conditions
trialset = [{subject(n).decision==0  & good}, ...
             {subject(n).decision==1  & good}];

% Precompute x-axis shift
x_shifted = (signal_bounds(1):signal_bounds(2)) - 6000;

% Initialize figure
figure; hold on;

for j = 1:length(trialset)
    % Get trial indices
    trialIDs = trialset{j}(1:end);

    % Extract high gamma matrix for all trials in one step
    high_gamma_mat = subject(n).electrode(e).trigger(t).high_gamma_mat(trialIDs, :);

    high_gamma_mat_baseline = subject(n).electrode(e).trigger(1).high_gamma_mat(trialIDs, :);
    % Compute baseline mean across trials
    baseline = mean(high_gamma_mat_baseline(:, baseline_bounds(1):baseline_bounds(2)), 2);

    % Normalize signal window (subtract baseline, divide by baseline mean)
    for_psth = (high_gamma_mat(:, signal_bounds(1):signal_bounds(2)) - baseline) ./ baseline;

    % Compute mean and standard error across trials (no smoothing)
    y = mean(for_psth, 1);
    dy = std(for_psth, [], 1) / sqrt(size(for_psth, 1));

    % Define fill colors based on trial index
    fill_color = [1 - (j - 1) / (length(trialset) - 0.9), (j - 1) / (length(trialset) - 0.9), 0.15];

    % Plot shaded error region (without smoothing)
    fill([x_shifted, fliplr(x_shifted)], [y - dy, fliplr(y + dy)], ...
         fill_color, 'linestyle', 'none', 'FaceAlpha', 0.6);
    
    % Plot mean signal (without smoothing)
    plot(x_shifted, y, 'LineWidth', 2, 'Color', [0.7*(1 - (j - 1) / (length(trialset) - 1)), ...
                                                 0.5*(j - 1) / (length(trialset) - 1), 0, 1]);
end

% Adjust plot appearance
set(gca, 'TickDir', 'out');
box off;
set(gcf, 'Color', 'w');
xlabel('Time from decision (milliseconds)', 'FontSize', 14);
ylabel('HFA', 'FontSize', 14);
xline(0);

%ylim([-0.05 0.25])


