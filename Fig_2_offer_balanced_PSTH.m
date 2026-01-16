%% Offer-weighted (by min #trials) PSTH: electrode-weighted (electrode = sample), need to run the Fig 2 code before running this
clear trialset
figure(1); clf

clust_indices = find(sorted_med_coor<5 & sorted_med_coor>-5); % MOS adjacent
%clust_indices = find(sorted_med_coor>(round(max(sorted_med_coor))-15)); % lat-most 1.5cm


t = 3;
num_conditions = 2; % 1=avoid, 2=approach

signal_bounds   = [4000 8000];
baseline_bounds = [5500 6000];
window = signal_bounds(2) - signal_bounds(1) + 1;

smoothing_kernel = 50;

% colors
red_line  = [230, 41, 41]/255;  red_fill  = [247, 189, 189]/255;
blue_line = [18, 6, 255]/255;   blue_fill = [180, 176, 255]/255;

% Store electrode-level PSTHs per condition (after weighted offer-averaging)
for_anova = cell(1,num_conditions);

for j = 1:num_conditions

    if j == 1
        line_color = red_line;
        fill_color = red_fill;
    else
        line_color = blue_line;
        fill_color = blue_fill;
    end

    % One row per electrode (electrode = sample)
    for_psth_elec = nan(length(clust_indices), window);%

    for i = 1:length(clust_indices)

        n = sorted_subject_ID(clust_indices(i));
        e = sorted_electrode_ID(clust_indices(i));

        good = subject(n).electrode(e).trigger(t).good_trials==1 & ...
               subject(n).electrode(e).trigger(1).good_trials==1;


        % trial types you want to include (your conflict criterion)
        conflicting_types = find(subject(n).conflict_trial_type>0);

        % weighted sums across trial types
        sum_trace_w = zeros(1, window);
        sum_w       = 0;

        for tt = conflicting_types(:)'

            mask_type  = (subject(n).trial_type_trial == tt) & good;
            mask_avoid = mask_type & (subject(n).decision == 0);
            mask_appr  = mask_type & (subject(n).decision == 1);

            nA = sum(mask_appr);
            nV = sum(mask_avoid);

            % need at least 2 of each (at least 2 approached; 2 avoided) to define mean traces for each
            % decision type
            if nA < 3 || nV <3
                continue
            end

            w = 1; % <-- equal weight per offer type

            % pick decision condition for this subplot
            if j == 1
                trialIDs = mask_avoid;
            else
                trialIDs = mask_appr;
            end

            % compute PSTH for this (electrode, trialtype, decision)
            temp = subject(n).electrode(e).trigger(t).high_gamma_mat( ...
                trialIDs, signal_bounds(1):signal_bounds(2));

            temp_baseline = subject(n).electrode(e).trigger(1).high_gamma_mat( ...
                trialIDs, baseline_bounds(1):baseline_bounds(2));

            baseline = mean(temp_baseline(:), 'omitnan');
            trace = (mean(temp,1,'omitnan') - baseline) ./ baseline;

            % accumulate weighted
            sum_trace_w = sum_trace_w + w * trace;
            sum_w       = sum_w + w;
        end

        % finalize electrode trace (weighted mean across trial types)
        if sum_w > 0
            for_psth_elec(i,:) = sum_trace_w ./ sum_w;
        end
    end

    % Keep only electrodes with data
    keep = ~isnan(for_psth_elec(:,1));
    for_psth = for_psth_elec(keep,:);

    % Mean + SEM across electrodes (electrode = sample)
    y  = mean(for_psth, 1, 'omitnan');
    dy = std(for_psth, 0, 1, 'omitnan') ./ sqrt(size(for_psth,1));

    for_anova{j} = for_psth;

    % Plot mean ± SEM
    time_ms = (signal_bounds(1):signal_bounds(2)) - 6000;
    fill([time_ms, fliplr(time_ms)], [y-dy, fliplr(y+dy)], ...
        fill_color, 'LineStyle','none', 'FaceAlpha', 0.6);
    hold on
    plot(time_ms, smooth(y, smoothing_kernel), 'LineWidth', 2, 'Color', line_color);
    hold on
end

% cosmetics
set(gca, 'TickDir', 'out');
box off; set(gcf, 'Color', 'w');
xlabel('time from decision (milliseconds)', 'FontSize', 20);
ylabel('\Delta HFA', 'Interpreter', 'tex', 'FontSize', 20);
xline(0,'k--','LineWidth',2)
xlim([-1800 500])

% electrode-level paired test in chosen window
idx = 1700:2000;
a = sum(for_anova{1}(:,idx), 2);  % avoid
b = sum(for_anova{2}(:,idx), 2);  % approach
keep = isfinite(a) & isfinite(b);
[h,p,ci,stats] = ttest(a(keep), b(keep))


%%

%% look at group PSTH
clear trialset
figure(1)
clust_indices = find(sorted_med_coor<5 & sorted_med_coor>-5); %MOS adjacent
%clust_indices = find(sorted_med_coor>(round(max(sorted_med_coor))-15) ); %lat-most 1.5cm
%clust_indices = find(sorted_lat_coor>0); %lat-most 1.5cm

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

% We will APPEND rows: each row = (electrode, trial_type)
for_psth_rows = [];   % [nRows x window]

% We will end with ONE row per electrode (after averaging across trial types)
for_psth_elec = nan(length(clust_indices), window);  % [nElectrodes x window]

for i = 1:length(clust_indices)

    n = sorted_subject_ID(clust_indices(i));
    e = sorted_electrode_ID(clust_indices(i));

    good = subject(n).electrode(e).trigger(t).good_trials==1 & ...
           subject(n).electrode(e).trigger(1).good_trials==1;

    conflicting_types = find(subject(n).conflict_trial_type>0);

    sum_trace = zeros(1, window);
    n_types   = 0;

    for tt = conflicting_types(:)'

        mask_type  = (subject(n).trial_type_trial == tt) & good;
        mask_avoid = mask_type & (subject(n).decision == 0);
        mask_appr  = mask_type & (subject(n).decision == 1);

        % require BOTH decisions exist for that type (so type is "eligible")
        if sum(mask_avoid)<3|| sum(mask_appr)<3
            continue
        end

        if j == 1
            trialIDs = mask_avoid;
        else
            trialIDs = mask_appr;
        end

        temp = subject(n).electrode(e).trigger(t).HG_agg_mat(trialIDs, signal_bounds(1):signal_bounds(2));
        temp_baseline = subject(n).electrode(e).trigger(1).HG_agg_mat(trialIDs, baseline_bounds(1):baseline_bounds(2));
        baseline = mean(temp_baseline(:), 'omitnan');

        trace = (mean(temp,1,'omitnan') - baseline) ./ baseline;

        % accumulate across trial types (equal weight per type)
        sum_trace = sum_trace + trace;
        n_types   = n_types + 1;
    end

    if n_types > 0
        for_psth_elec(i,:) = sum_trace ./ n_types;  % average across types
    end
end

% =========================
% Subject-equal weighting
% =========================

% Identify which subject each electrode row came from
subj_of_row = sorted_subject_ID(clust_indices);   % length = nElectrodes
subj_of_row = subj_of_row(:);

% Keep only electrodes with data
keep = ~isnan(for_psth_elec(:,1));
E = for_psth_elec(keep,:);
subj_keep = subj_of_row(keep);

uSubj = unique(subj_keep);

% One row per subject: average across that subject's electrodes
subj_psth = nan(numel(uSubj), window);

for s = 1:numel(uSubj)
    idx = subj_keep == uSubj(s);
    subj_psth(s,:) = mean(E(idx,:), 1, 'omitnan');
end

% Plot group mean/SEM across subjects (each subject = 1 vote)
y  = mean(subj_psth, 1, 'omitnan');
dy = std(subj_psth, 0, 1, 'omitnan') ./ sqrt(size(subj_psth,1));

% Store for later stats
for_anova{j} = subj_psth;   % now rows = subjects, NOT electrodes




        fill([(signal_bounds(1):signal_bounds(2))-6000, fliplr((signal_bounds(1):signal_bounds(2))-6000)], [(y - dy)  fliplr((y + dy))],fill_color, 'linestyle', 'none','FaceAlpha', 0.6);
        hold on
        plot((signal_bounds(1):signal_bounds(2))'-6000, smooth(y,smoothing_kernel), 'Linewidth',2,'Color',line_color)
        hold on


end

% Adjust plot appearance
set(gca, 'TickDir', 'out');
% Get the current x-axis limits
x_limits = xlim;


% Get the current y-axis limits
y_limits = ylim;

% Define tick marks, ensuring that 0 is included
y_ticks = linspace(y_limits(1), y_limits(2), 2);  % Create 5 evenly spaced ticks
y_ticks = unique([y_ticks, 0]);  % Ensure 0 is included

% Set the y-axis ticks
set(gca, 'YTick', y_ticks);


box off;
set(gcf, 'Color', 'w');
xlabel('time from decision (milliseconds)', 'FontSize', 20);
ylabel('\Delta HFA', 'Interpreter', 'tex', 'FontSize', 20);

xline(0,'k--','Linewidth',2)
xlim([-2000 500])

[h p c stats]=ttest([sum(for_anova{1}(:,1000:2000),2)'],[sum(for_anova{2}(:,1000:2000),2)'])

%% Offer-balanced PSTH: approach vs avoid
clear trialset
figure(1); clf
clust_indices = find(sorted_med_coor<5 & sorted_med_coor>-5); % MOS adjacent
%clust_indices = find(sorted_med_coor>(round(max(sorted_med_coor))-15) ); %lat-most 1.5cm


t = 3; num_conditions = 2;  % 1=avoid, 2=approach

signal_bounds   = [3000 8000];
baseline_bounds = [5500 6000];
window = signal_bounds(2) - signal_bounds(1) + 1;


% colors
red_line  = [230, 41, 41]/255;  red_fill  = [247, 189, 189]/255;
blue_line = [18, 6, 255]/255;   blue_fill = [180, 176, 255]/255;

% Map each clust row -> subject id
subj_of_row = sorted_subject_ID(clust_indices(:));
uSubj = unique(subj_of_row);

% Store subject-level PSTHs per condition
subj_psth = cell(1,num_conditions);
for j = 1:num_conditions
    subj_psth{j} = nan(numel(uSubj), window);
end

% ---------------------------
% Build subject-level PSTH
% ---------------------------
for s = 1:numel(uSubj)
    subj = uSubj(s);

    % electrode rows belonging to this subject
    rows_s = find(subj_of_row == subj);

    % these will accumulate electrode-level offer-balanced traces
    elec_traces = cell(1,num_conditions);
    for j = 1:num_conditions
        elec_traces{j} = [];  % will become [nElecUsed x window]
    end

    for rr = 1:numel(rows_s)
        idx = rows_s(rr);
        n = sorted_subject_ID(clust_indices(idx));
        e = sorted_electrode_ID(clust_indices(idx));

        good = subject(n).electrode(e).trigger(t).good_trials==1 & ...
               subject(n).electrode(e).trigger(1).good_trials==1;

        % define "offer type" as (reward_offer, punishment_offer) pair
        % using your existing trial-level vectors:
        R = subject(n).trial_type_trial(:);
        P = subject(n).punish_trial(:);

        % only consider trials that are "conflictful" if you want:
        % (if you want ALL offers, comment this line out)
        ok_conflict = subject(n).conflict_trial(:)  > 0;

        % include only valid trials
        base_mask = good(:)  & ok_conflict & subject(n).decision(:) < 2;

        offers = unique([R(base_mask) P(base_mask)], 'rows');

        % For THIS electrode: accumulate per-offer traces (equal weight per offer)
        sum_trace = cell(1,num_conditions);
        n_offers  = 0;
        for j = 1:num_conditions
            sum_trace{j} = zeros(1, window);
        end
        total_weight=0;

        for k = 1:size(offers,1)
            r0 = offers(k,1);
            p0 = offers(k,2);

            mask_offer = base_mask & (R == r0) & (P == p0);

            mask_avoid = mask_offer & (subject(n).decision(:) == 0);
            mask_appr  = mask_offer & (subject(n).decision(:) == 1);


            % weight for this offer = min(#avoid, #approach)
            w = min(sum(mask_avoid), sum(mask_appr));
            if w <= 1, continue; end

            w = 1;

            % -------- avoid trace (j=1) --------
            trialIDs = mask_avoid;
            temp = subject(n).electrode(e).trigger(t).high_gamma_mat(trialIDs, signal_bounds(1):signal_bounds(2));
            temp_baseline = subject(n).electrode(e).trigger(1).high_gamma_mat(trialIDs, baseline_bounds(1):baseline_bounds(2));
            baseline = mean(temp_baseline(:), 'omitnan');
            trace_avoid = (mean(temp,1,'omitnan') - baseline) ./ baseline;

            % -------- approach trace (j=2) --------
            trialIDs = mask_appr;
            temp = subject(n).electrode(e).trigger(t).high_gamma_mat(trialIDs, signal_bounds(1):signal_bounds(2));
            temp_baseline = subject(n).electrode(e).trigger(1).high_gamma_mat(trialIDs, baseline_bounds(1):baseline_bounds(2));
            baseline = mean(temp_baseline(:), 'omitnan');
            trace_appr = (mean(temp,1,'omitnan') - baseline) ./ baseline;

            % accumulate weighted sums (same w for both)
            sum_trace{1} = sum_trace{1} + w * trace_avoid;
            sum_trace{2} = sum_trace{2} + w * trace_appr;
            total_weight = total_weight + w;

            n_offers = n_offers + 1;

        end

        % finalize electrode-level offer-balanced trace (one per decision)
        if total_weight > 0
            elec_traces{1}(end+1,:) = sum_trace{1} ./ total_weight;
            elec_traces{2}(end+1,:) = sum_trace{2} ./ total_weight;
        end

    end

    % Now average electrodes within subject (subject-equal weighting)
    for j = 1:num_conditions
        if ~isempty(elec_traces{j})
            subj_psth{j}(s,:) = mean(elec_traces{j}, 1, 'omitnan');
        end
    end
end

% ---------------------------
% Plot group mean ± SEM across subjects
% ---------------------------
time_ms = (signal_bounds(1):signal_bounds(2)) - 6000;  % your convention

for j = 1:num_conditions
    if j == 1
        line_color = red_line;  fill_color = red_fill;
    else
        line_color = blue_line; fill_color = blue_fill;
    end

    X = subj_psth{j};
    keepS = ~isnan(X(:,1));
    X = X(keepS,:);

    y  = mean(X, 1, 'omitnan');
    dy = std(X, 0, 1, 'omitnan') ./ sqrt(size(X,1));

    fill([time_ms, fliplr(time_ms)], [y-dy, fliplr(y+dy)], ...
        fill_color, 'LineStyle','none', 'FaceAlpha', 0.6);
    hold on
    plot(time_ms, y, 'LineWidth', 2, 'Color', line_color);
    hold on
end

box off; set(gca,'TickDir','out'); set(gcf,'Color','w');
xlabel('time from decision (ms)', 'FontSize', 20);
ylabel('\Delta HFA', 'FontSize', 20);
xline(0,'k--','LineWidth',2);
xlim([-1500 500]);
%%
% ---------------------------
% Plot group mean ± SEM across subjects (smooth mean trace for display)
% ---------------------------
time_ms = (signal_bounds(1):signal_bounds(2)) - 6000;  % your convention
smooth_ms = 50;   % 50 ms smoothing for display

for j = 1:num_conditions
    if j == 1
        line_color = red_line;  fill_color = red_fill;
    else
        line_color = blue_line; fill_color = blue_fill;
    end

    X = subj_psth{j};
    keepS = ~isnan(X(:,1));
    X = X(keepS,:);

    y  = mean(X, 1, 'omitnan');
    dy = std(X, 0, 1, 'omitnan') ./ sqrt(size(X,1));

    % fill uses unsmoothed mean/SEM (recommended)
    fill([time_ms, fliplr(time_ms)], [y-dy, fliplr(y+dy)], ...
        fill_color, 'LineStyle','none', 'FaceAlpha', 0.6);
    hold on

    % smooth only the plotted line (moving average over 50 ms)
    y_smooth = smoothdata(y, 'movmean', smooth_ms);

    plot(time_ms, y_smooth, 'LineWidth', 2, 'Color', line_color);
    hold on
end

box off; set(gca,'TickDir','out'); set(gcf,'Color','w');
xlabel('time from decision (ms)', 'FontSize', 20);
ylabel('\Delta HFA', 'FontSize', 20);
xline(0,'k--','LineWidth',2);
xlim([-2000 500]);

%% anova
%for the data aligned to decision
spacing = 50; start_sliding = 500; length_epoch = 500; 
last_sliding = 3000 - length_epoch; num_epochs = (last_sliding - start_sliding) / spacing + 1;

bins = zeros(num_epochs, 2);
bins(1, :) = [start_sliding (start_sliding + length_epoch - 1)];

for i = 2:num_epochs
    bins(i, :) = [(bins(i-1, 1) + spacing) (bins(i-1, 2) + spacing)];
end

t_stat = []; ps=[];psr=[];
for i = 1:num_epochs
    [h p c stats]=ttest([sum(subj_psth{1}(:,bins(i,1):bins(i,2)),2)],[sum(subj_psth{2}(:,bins(i,1):bins(i,2)),2)])
    t_stat = [t_stat stats.tstat];
    ps = [ps p]
    p_signrank = signrank([sum(subj_psth{1}(:,bins(i,1):bins(i,2)),2)'],[sum(subj_psth{2}(:,bins(i,1):bins(i,2)),2)']);     % paired, nonparametric
    psr = [psr p_signrank]
end
%%
ps
[h p c stats]=ttest([sum(subj_psth{1}(:,2400:2800),2)'],[sum(subj_psth{2}(:,2400:3000),2)'])
p_signrank = signrank([sum(subj_psth{1}(:,2500:3000),2)'],[sum(subj_psth{2}(:,2500:3000),2)'])     % paired, nonparametric
