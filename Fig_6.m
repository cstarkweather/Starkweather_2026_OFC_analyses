%% Figure 3G: 10ms bin occupancy, SE patches visible, stats on 100ms windows
figure(2); clf

% ----------------- params -----------------
condition      = 2;     % 1 = approach, 2 = avoid
spacing        = 5;     % ms step between bins
start_sliding  = -1500;
length_epoch   = 10;    % BIN size for occupancy (10 ms)
last_sliding   = 0 ;
num_bins       = (last_sliding - start_sliding)/spacing + 1;
smooth_ms   = 10;                          % smoothing window in ms (edit)
alpha          = 0.05;

% stats window settings (preferential occupancy over 100 ms)
win_stat_ms    = 200;   % window length for stats
step_stat_ms   = 50;     % step for stats windows (often same as spacing)

red_line = [230, 41, 41]/255;
red_fill = [247, 189, 189]/255;

blue_line = [18, 6, 255]/255;
blue_fill = [180, 176, 255]/255;

% ----------------- build 10ms bins -----------------
bins = zeros(num_bins,2);
bins(1,:) = [start_sliding, start_sliding + length_epoch - 1];
for i = 2:num_bins
    bins(i,:) = bins(i-1,:) + spacing;
end

% x-axis you want to plot (your convention)
time_vec = bins(:,1) + length_epoch;  % shifts by window length

% ----------------- storage across subjects -----------------
p1_subj = nan(num_bins,6);
p2_subj = nan(num_bins,6);

total_p1 = zeros(num_bins,1);
total_p2 = zeros(num_bins,1);
total_valid = zeros(num_bins,1);

% ----------------- compute occupancy per subject -----------------
for n = 1:6
    load(sprintf('state_metadata_%d.mat', n)); % entry_times, departure_times, trial_identity, state_identity

    % most_conflicting trial type
    most_conflicting = find(subject(n).conflict_trial_type == max(subject(n).conflict_trial_type));
    if numel(most_conflicting) > 1
        winner = find(abs(subject(n).value_trial_type(most_conflicting)) == ...
                      min(abs(subject(n).value_trial_type(most_conflicting))), 1, 'first');
        most_conflicting = most_conflicting(winner);
    end

    % ONE condition logic
    switch condition
        case 1
            filter = find(subject(n).electrode(1).trigger(3).good_trials == 1 & ...
                          subject(n).decision == 1 & ...
                          subject(n).trial_type_trial == most_conflicting);    
                          %subject(n).reward_trial == 7 & subject(n).punish_trial == 0);
                                     
        case 2
            filter = find(subject(n).electrode(1).trigger(3).good_trials == 1 & ...
                          subject(n).decision == 0 & ...
                          subject(n).trial_type_trial == most_conflicting);   
                          %subject(n).trial_type_trial == length(unique(subject(n).trial_type_trial)));
                               
    end
    if isempty(filter), continue; end

    p_state = zeros(num_bins,2);
    valid   = zeros(num_bins,1);

    for ii = 1:numel(filter)
        tr = filter(ii);

        entry_t = entry_times(trial_identity == tr);
        dep_t   = departure_times(trial_identity == tr);
        st_id   = state_identity(trial_identity == tr);

        rt_tr = subject(n).rt(tr);

        for j = 1:num_bins
            % RT gating
            if rt_tr <= -(bins(j,2))
                continue;
            end
            valid(j) = valid(j) + 1;

            bin_ind = find(entry_t < bins(j,1) & (dep_t > bins(j,2) | isnan(dep_t)));
            for a = 1:numel(bin_ind)
                if st_id(bin_ind(a)) == 1
                    p_state(j,1) = p_state(j,1) + 1;
                else
                    p_state(j,2) = p_state(j,2) + 1;
                end
            end
        end
    end

    p1 = p_state(:,1)./valid;
    p2 = p_state(:,2)./valid;
    p1(valid==0) = NaN;  p2(valid==0) = NaN;

    p1_subj(:,n) = p1;
    p2_subj(:,n) = p2;

    total_p1    = total_p1 + p1 .* valid;
    total_p2    = total_p2 + p2 .* valid;
    total_valid = total_valid + valid;
end

% ----------------- SMOOTHING (subject-level, NaN-safe) -----------------
smooth_bins = max(1, round(smooth_ms/spacing));

p1_subj_sm = p1_subj;
p2_subj_sm = p2_subj;

for n = 1:6
    p1_subj_sm(:,n) = smooth_nan_segments(p1_subj(:,n), smooth_bins);
    p2_subj_sm(:,n) = smooth_nan_segments(p2_subj(:,n), smooth_bins);
end

% ----------------- group mean + SE (binwise) on SMOOTHED traces -----------------
avg1 = nanmean(p1_subj_sm, 2);
avg2 = nanmean(p2_subj_sm, 2);

se1  = nanstd(p1_subj_sm, 0, 2) ./ sqrt(6);
se2  = nanstd(p2_subj_sm, 0, 2) ./ sqrt(6);

% If you still want the weighted mean by valid bins (optional), keep your old avg calc.
% In that case, only use smoothing for visualization (avg/se), and keep stats on raw.




% ----------------- stats on 100ms windows built from 10ms-bin series -----------------
% build 100ms windows over the same time axis (use bin START times for inclusion)
t_bin_start = bins(:,1);                 % bin start times (ms)
t_bin_end   = bins(:,2);                 % bin end times (ms)

stat_starts = (min(t_bin_start):step_stat_ms:(0 - win_stat_ms))';  % e.g., -1500 to -100
stat_ends   = stat_starts + win_stat_ms - 1;
num_w = numel(stat_starts);

p_right = nan(num_w,1); % state1 > state2
p_left  = nan(num_w,1); % state2 > state1

for w = 1:num_w
    w0 = stat_starts(w);
    w1 = stat_ends(w);

    % bins fully contained in this 100ms window
    inwin = (t_bin_start >= w0) & (t_bin_end <= w1);

    % per-subject window-averaged occupancy
    p1w = nanmean(p1_subj(inwin,:), 1)';  % [6x1]
    p2w = nanmean(p2_subj(inwin,:), 1)';  % [6x1]
    d   = p1w - p2w;

    d = d(~isnan(d));
    if numel(d) < 3, continue; end

    [~,p_right(w)] = ttest(d,0,'Tail','right');
    [~,p_left(w)]  = ttest(d,0,'Tail','left');
end

sig12 = p_right < alpha;
sig21 = p_left  < alpha;

% convert significant windows to contiguous ranges (in ms)
ranges12 = local_sig_ranges(sig12, [stat_starts stat_ends]);
ranges21 = local_sig_ranges(sig21, [stat_starts stat_ends]);

% ----------------- plotting -----------------
hold on

% 1) significance shading FIRST (background)
yl = [0 0.5];
ylim(yl);
xlim([-1000 0]);

t_shift = 0;

for k = 1:numel(ranges12)
    t0 = ranges12{k}(1) + t_shift;
    t1 = ranges12{k}(2) + t_shift;
    h = patch([t0 t1 t1 t0], [yl(1) yl(1) yl(2) yl(2)], blue_fill, ...
              'FaceAlpha', 0.10, 'EdgeColor','none');
    set(h,'Tag','sigPatch');
end

for k = 1:numel(ranges21)
    t0 = ranges21{k}(1) + t_shift;
    t1 = ranges21{k}(2) + t_shift;
    h = patch([t0 t1 t1 t0], [yl(1) yl(1) yl(2) yl(2)], red_fill, ...
              'FaceAlpha', 0.10, 'EdgeColor','none');
    set(h,'Tag','sigPatch');
end

% 2) SE patches (draw only where non-NaN so they actually show)
hSe1 = fill_nan_gaps(time_vec, avg1 - se1, avg1 + se1, blue_fill, 0.30);
hSe2 = fill_nan_gaps(time_vec, avg2 - se2, avg2 + se2, red_fill,  0.30);
set(hSe1,'Tag','seFill'); set(hSe2,'Tag','seFill');

% 3) mean lines on top
l1 = plot(time_vec, avg1, 'Color', blue_line, 'LineWidth', 2);
l2 = plot(time_vec, avg2, 'Color', red_line,  'LineWidth', 2);

% enforce stacking
uistack(findobj(gca,'Tag','sigPatch'),'bottom');
uistack(findobj(gca,'Tag','seFill'),'top');
uistack([l1 l2],'top');

% style
set(gca,'TickDir','out'); box off;
set(gcf,'Color','w');
xlabel('Time (ms) relative to decision','FontSize',20);
ylabel('p(state)','FontSize',24);

% print ranges
fprintf('Significant 100ms windows where p(state1) > p(state2) (uncorrected, alpha=%.3f):\n', alpha);
if isempty(ranges12), disp('  none'); else
    for k = 1:numel(ranges12)
        fprintf('  %4d to %4d ms\n', ranges12{k}(1)+t_shift, ranges12{k}(2)+t_shift);
    end
end

fprintf('Significant 100ms windows where p(state2) > p(state1) (uncorrected, alpha=%.3f):\n', alpha);
if isempty(ranges21), disp('  none'); else
    for k = 1:numel(ranges21)
        fprintf('  %4d to %4d ms\n', ranges21{k}(1)+t_shift, ranges21{k}(2)+t_shift);
    end
end


%% Count state ENTRIES in last 400ms (state1 vs state2), approach vs avoid
figure; clf

w0 = -400;     % last 400 ms pre-decision
%t0 = -last_window_ms;
t1 = 0;

counts_appr = nan(6,2);   % [subj x (state1,state2)]  (RAW counts)
counts_avoid = nan(6,2);

ntr_appr  = nan(6,1);     % number of valid trials per subject (approach)
ntr_avoid = nan(6,1);     % number of valid trials per subject (avoid)

for n = [1 4 5 6]
    Y = [];           % 1 = state1, 0 = state2
Decision = [];    % 1 = approach, 0 = avoid
Subject = [];

    load(sprintf('state_metadata_%d.mat', n)); % entry_times, departure_times, trial_identity, state_identity

    % highest conflict trial type (your logic)
    most_conflicting = find(subject(n).conflict_trial_type == max(subject(n).conflict_trial_type));
    if numel(most_conflicting) > 1
        winner = find(abs(subject(n).value_trial_type(most_conflicting)) == ...
                      min(abs(subject(n).value_trial_type(most_conflicting))), 1, 'first');
        most_conflicting = most_conflicting(winner);
    end

    for condition = 1:2
        switch condition
            case 1 % approach
                filter = find(subject(n).electrode(1).trigger(3).good_trials == 1 & ...
                              subject(n).decision == 1 & ...
                              subject(n).trial_type_trial == most_conflicting);
            case 2 % avoid
                filter = find(subject(n).electrode(1).trigger(3).good_trials == 1 & ...
                              subject(n).decision == 0 & ...
                              subject(n).trial_type_trial == most_conflicting);
        end
        if isempty(filter), continue; end

        c1 = 0; c2 = 0;
        ntr = 0;

        for ii = 1:numel(filter)
            tr = filter(ii);

            ntr = ntr + 1;

            % states for this trial
            e_t  = entry_times(trial_identity == tr);
            d_t  = departure_times(trial_identity == tr);
            s_id = state_identity(trial_identity == tr);

            % treat open-ended departures as +inf
            d_t(isnan(d_t)) = inf;

            % (A) count state active at w0
            active = find(e_t <= w0 & d_t > w0);
            if ~isempty(active)
                k0 = active(end);
                if s_id(k0) == 1
                    c1 = c1 + 1;
                else
                    c2 = c2 + 1;
                end
            end

            % (B) count new entries after w0
            inwin_entries = (e_t > w0) & (e_t <= t1);
            c1 = c1 + sum(inwin_entries & (s_id == 1));
            c2 = c2 + sum(inwin_entries & (s_id ~= 1));  % assumes only states 1/2
        end

        % store raw counts + ntr
        if condition == 1
            counts_appr(n,:) = [c1 c2];
            ntr_appr(n)      = ntr;
        else
            counts_avoid(n,:) = [c1 c2];
            ntr_avoid(n)      = ntr;
        end
    end
end

% Normalize by # trials (entries per trial)
norm_appr  = counts_appr ./ ntr_appr;    % implicit expansion: [6x2] ./ [6x1]
norm_avoid = counts_avoid ./ ntr_avoid;



% Plot: each subject is a line from state1 -> state2 (normalized)
subplot(1,2,1); hold on
for n = 1:6
    if any(isnan(norm_appr(n,:))) || isnan(ntr_appr(n)) || ntr_appr(n)==0, continue; end
    plot([1 2], norm_appr(n,:), '-', 'LineWidth', 1.5,'Color',[(n-1)/6 (n-1)/6 (n-1)/6]);
end
set(gca,'XTick',[1 2],'XTickLabel',{'State 1 entries/trial','State 2 entries/trial'},'TickDir','out');
title('Approach (highest conflict): last 400 ms'); ylabel('Entries per trial'); box off

subplot(1,2,2); hold on
for n = 1:6
    if any(isnan(norm_avoid(n,:))) || isnan(ntr_avoid(n)) || ntr_avoid(n)==0, continue; end
    plot([1 2], norm_avoid(n,:), '-', 'LineWidth', 1.5,'Color',[(n-1)/6 (n-1)/6 (n-1)/6]);
end
set(gca,'XTick',[1 2],'XTickLabel',{'State 1 entries/trial','State 2 entries/trial'},'TickDir','out');
title('Avoid (highest conflict): last 400 ms'); ylabel('Entries per trial'); box off

set(gcf,'Color','w');


%%
%% ============================================================
% TRIAL-LEVEL ANALYSIS:
% Does decision affect state-entry preference on a PER-TRIAL basis? For the
% highest level conflict trial only; restricting to subject where >2 of
% each choice in each trial type
% Preference per trial = (#State1 entries) - (#State2 entries)
% Model: TrialPref ~ Decision + (1|Subject)
%% ============================================================

TrialPref = [];
DecisionT = [];
SubjectT  = [];
TrialT    = [];

w0 = -400;
t1 = 0;

for n = [1 4 5 6]
    load(sprintf('state_metadata_%d.mat', n)); % entry_times, departure_times, trial_identity, state_identity

    % highest conflict trial type (same logic)
    most_conflicting = find(subject(n).conflict_trial_type == max(subject(n).conflict_trial_type));
    if numel(most_conflicting) > 1
        winner = find(abs(subject(n).value_trial_type(most_conflicting)) == ...
                      min(abs(subject(n).value_trial_type(most_conflicting))), 1, 'first');
        most_conflicting = most_conflicting(winner);
    end

    % valid trials for this subject
    valid_trials = find(subject(n).electrode(1).trigger(3).good_trials == 1 & ...
                        subject(n).trial_type_trial == most_conflicting & ...
                        subject(n).decision < 2);

    for ii = 1:numel(valid_trials)
        tr = valid_trials(ii);

        % decision for this trial
        decVal = subject(n).decision(tr);   % 1=approach, 0=avoid

        % state segments for this trial
        e_t  = entry_times(trial_identity == tr);
        d_t  = departure_times(trial_identity == tr);
        s_id = state_identity(trial_identity == tr);

        d_t(isnan(d_t)) = inf;

        c1 = 0; c2 = 0;

        % (A) state active at w0
        active = find(e_t <= w0 & d_t > w0);
        if ~isempty(active)
            k0 = active(end);
            if s_id(k0) == 1
                c1 = c1 + 1;
            else
                c2 = c2 + 1;
            end
        end

        % (B) new entries after w0
        inwin_entries = (e_t > w0) & (e_t <= t1);
        c1 = c1 + sum(inwin_entries & (s_id == 1));
        c2 = c2 + sum(inwin_entries & (s_id ~= 1));

        % trial-level preference
        TrialPref(end+1,1) = c1 - c2;
        DecisionT(end+1,1) = decVal;
        SubjectT(end+1,1)  = n;
        TrialT(end+1,1)    = tr;
    end
end

fprintf('\nTrial-level dataset:\n');
fprintf('  N trials: %d\n', numel(TrialPref));
fprintf('  Mean TrialPref: %.3f\n', mean(TrialPref));

% ------------------------------------------------------------
% Mixed-effects model (Gaussian; counts can be negative)
% ------------------------------------------------------------
tblT = table(TrialPref, DecisionT, categorical(SubjectT), TrialT, ...
    'VariableNames', {'TrialPref','Decision','Subject','Trial'});

lme = fitlme(tblT, ...
    'TrialPref ~ Decision + (1|Subject)');

disp('--- LME: TrialPref ~ Decision + (1|Subject) ---');
disp(lme);

%%
%% Figure 3G (single subject): 10ms bin occupancy, trial-SE patches, stats on 200ms windows
figure(2); clf

% ----------------- params -----------------
n0             = 4;     % <-- choose subject here
condition      = 1;     % 1 = approach, 2 = avoid
spacing        = 5;     % ms step between bins
start_sliding  = -1500;
length_epoch   = 10;    % BIN size for occupancy (10 ms)
last_sliding   = 0;
num_bins       = (last_sliding - start_sliding)/spacing + 1;

smooth_ms      = 10;    % smoothing window (ms) for visualization only
alpha          = 0.05;

% stats window settings (preferential occupancy over win_stat_ms)
win_stat_ms    = 200;
step_stat_ms   = 50;

red_line  = [230, 41, 41]/255;
red_fill  = [247, 189, 189]/255;
blue_line = [18, 6, 255]/255;
blue_fill = [180, 176, 255]/255;

% ----------------- build 10ms bins -----------------
bins = zeros(num_bins,2);
bins(1,:) = [start_sliding, start_sliding + length_epoch - 1];
for i = 2:num_bins
    bins(i,:) = bins(i-1,:) + spacing;
end

time_vec = bins(:,1) + length_epoch;  % your convention

% ----------------- load state metadata for subject -----------------
load(sprintf('state_metadata_%d.mat', n0)); % entry_times, departure_times, trial_identity, state_identity

% ----------------- choose most conflicting trial type (same as yours) -----------------
most_conflicting = find(subject(n0).conflict_trial_type == max(subject(n0).conflict_trial_type));
if numel(most_conflicting) > 1
    winner = find(abs(subject(n0).value_trial_type(most_conflicting)) == ...
                  min(abs(subject(n0).value_trial_type(most_conflicting))), 1, 'first');
    most_conflicting = most_conflicting(winner);
end

% ----------------- trial filter for this subject -----------------
switch condition
    case 1
        filter = find(subject(n0).electrode(1).trigger(3).good_trials == 1 & ...
                      subject(n0).decision == 1 & ...
                      subject(n0).trial_type_trial == most_conflicting);
    case 2
        filter = find(subject(n0).electrode(1).trigger(3).good_trials == 1 & ...
                      subject(n0).decision == 0 & ...
                      subject(n0).trial_type_trial == most_conflicting);
end
if isempty(filter)
    error('No trials matched filter for subject %d.', n0);
end

nTr = numel(filter);

% ----------------- trial-by-trial occupancy matrices -----------------
% p1_tr(:,t) = fraction of time (0/1 here) bin is in state1 for trial t
p1_tr = nan(num_bins, nTr);
p2_tr = nan(num_bins, nTr);

for ii = 1:nTr
    tr = filter(ii);

    entry_t = entry_times(trial_identity == tr);
    dep_t   = departure_times(trial_identity == tr);
    st_id   = state_identity(trial_identity == tr);

    rt_tr = subject(n0).rt(tr);

    for j = 1:num_bins
        % RT gating (same rule you used)
        if rt_tr <= -(bins(j,2))
            continue;
        end

        % states overlapping this bin
        bin_ind = find(entry_t < bins(j,1) & (dep_t > bins(j,2) | isnan(dep_t)));

        % if multiple segments overlap, mark as occupied if ANY overlap is that state
        % (this mimics your counting; if you want "fraction of overlap", that’s a different metric)
        has1 = any(st_id(bin_ind) == 1);
        has2 = any(st_id(bin_ind) ~= 1);

        p1_tr(j,ii) = double(has1);
        p2_tr(j,ii) = double(has2);
    end
end

% ----------------- optional smoothing (within-trial mean trace only) -----------------
smooth_bins = max(1, round(smooth_ms/spacing));

avg1_raw = nanmean(p1_tr, 2);
avg2_raw = nanmean(p2_tr, 2);

avg1 = smooth_nan_segments(avg1_raw, smooth_bins);
avg2 = smooth_nan_segments(avg2_raw, smooth_bins);

% SE across trials (compute on UNSMOOTHED trial series, then smooth SE only for display)
se1_raw = nanstd(p1_tr, 0, 2) ./ sqrt(sum(isfinite(p1_tr),2));
se2_raw = nanstd(p2_tr, 0, 2) ./ sqrt(sum(isfinite(p2_tr),2));

se1 = smooth_nan_segments(se1_raw, smooth_bins);
se2 = smooth_nan_segments(se2_raw, smooth_bins);

% ----------------- stats: within-subject across trials in sliding windows -----------------
t_bin_start = bins(:,1);
t_bin_end   = bins(:,2);

stat_starts = (min(t_bin_start):step_stat_ms:(0 - win_stat_ms))';
stat_ends   = stat_starts + win_stat_ms - 1;
num_w = numel(stat_starts);

p_right = nan(num_w,1); % state1 > state2
p_left  = nan(num_w,1); % state2 > state1

for w = 1:num_w
    w0 = stat_starts(w);
    w1 = stat_ends(w);

    inwin = (t_bin_start >= w0) & (t_bin_end <= w1);

    % per-trial window-averaged occupancy
    p1w = nanmean(p1_tr(inwin,:), 1)';  % [nTr x 1]
    p2w = nanmean(p2_tr(inwin,:), 1)';  % [nTr x 1]
    d   = p1w - p2w;

    d = d(isfinite(d));
    if numel(d) < 5, continue; end  % require enough trials

    [~,p_right(w)] = ttest(d, 0, 'Tail','right');
    [~,p_left(w)]  = ttest(d, 0, 'Tail','left');
end

sig12 = p_right < alpha;
sig21 = p_left  < alpha;

ranges12 = local_sig_ranges(sig12, [stat_starts stat_ends]);
ranges21 = local_sig_ranges(sig21, [stat_starts stat_ends]);

% ----------------- plotting -----------------
hold on
%yl = [0 0.5];
%ylim(yl);
xlim([-1000 0]);

% significance patches first
for k = 1:numel(ranges12)
    t0 = ranges12{k}(1);
    t1 = ranges12{k}(2);
    h = patch([t0 t1 t1 t0], [yl(1) yl(1) yl(2) yl(2)], blue_fill, ...
              'FaceAlpha', 0.10, 'EdgeColor','none');
    set(h,'Tag','sigPatch');
end

for k = 1:numel(ranges21)
    t0 = ranges21{k}(1);
    t1 = ranges21{k}(2);
    h = patch([t0 t1 t1 t0], [yl(1) yl(1) yl(2) yl(2)], red_fill, ...
              'FaceAlpha', 0.10, 'EdgeColor','none');
    set(h,'Tag','sigPatch');
end

% SE fills (only where non-NaN)
hSe1 = fill_nan_gaps(time_vec, avg1 - se1, avg1 + se1, blue_fill, 0.30);
hSe2 = fill_nan_gaps(time_vec, avg2 - se2, avg2 + se2, red_fill,  0.30);
set(hSe1,'Tag','seFill'); set(hSe2,'Tag','seFill');

% mean lines
l1 = plot(time_vec, avg1, 'Color', blue_line, 'LineWidth', 2);
l2 = plot(time_vec, avg2, 'Color', red_line,  'LineWidth', 2);

uistack(findobj(gca,'Tag','sigPatch'),'bottom');
uistack(findobj(gca,'Tag','seFill'),'top');
uistack([l1 l2],'top');

set(gca,'TickDir','out'); box off;
set(gcf,'Color','w');
xlabel('Time (ms) relative to decision','FontSize',20);
ylabel('p(state)','FontSize',24);
title(sprintf('Subject %d | condition %d | N trials=%d', n0, condition, nTr));

% print ranges
fprintf('Subject %d: significant %dms windows where p(state1) > p(state2) (alpha=%.3f):\n', ...
    n0, win_stat_ms, alpha);
if isempty(ranges12), disp('  none'); else
    for k = 1:numel(ranges12)
        fprintf('  %4d to %4d ms\n', ranges12{k}(1), ranges12{k}(2));
    end
end

fprintf('Subject %d: significant %dms windows where p(state2) > p(state1) (alpha=%.3f):\n', ...
    n0, win_stat_ms, alpha);
if isempty(ranges21), disp('  none'); else
    for k = 1:numel(ranges21)
        fprintf('  %4d to %4d ms\n', ranges21{k}(1), ranges21{k}(2));
    end
end


%% ----------------- helpers -----------------
function ranges = local_sig_ranges(sigmask, winmat)
% sigmask: [num_w x 1] logical
% winmat:  [num_w x 2] where each row is [start end]
    sig_idx = find(sigmask);
    ranges = {};
    if isempty(sig_idx), return; end

    d = diff(sig_idx);
    run_starts = [1; find(d > 1) + 1];
    run_ends   = [find(d > 1); numel(sig_idx)];

    for kk = 1:numel(run_starts)
        w1 = sig_idx(run_starts(kk));
        w2 = sig_idx(run_ends(kk));
        ranges{kk,1} = [winmat(w1,1), winmat(w2,2)];
    end
end

function h = fill_nan_gaps(x, ylo, yhi, faceColor, faceAlpha)
% Draws fill patches only over contiguous segments where both bounds are finite.
    good = isfinite(x) & isfinite(ylo) & isfinite(yhi);
    h = gobjects(0);

    if ~any(good), return; end
    idx = find(good);
    breaks = [1; find(diff(idx) > 1) + 1; numel(idx)+1];

    for b = 1:numel(breaks)-1
        seg = idx(breaks(b):breaks(b+1)-1);
        xs = x(seg);
        yl = ylo(seg);
        yh = yhi(seg);

        hh = fill([xs; flipud(xs)], [yl; flipud(yh)], faceColor, ...
                  'FaceAlpha', faceAlpha, 'EdgeColor','none');
        h(end+1) = hh; %#ok<AGROW>
    end
end

%% ===================== helper: smooth contiguous non-NaN segments =====================
function ysm = smooth_nan_segments(y, win)
    ysm = y;
    good = isfinite(y);
    if ~any(good), return; end

    idx = find(good);
    breaks = [1; find(diff(idx) > 1) + 1; numel(idx)+1];

    for b = 1:numel(breaks)-1
        seg = idx(breaks(b):breaks(b+1)-1);
        if numel(seg) < 3
            ysm(seg) = y(seg);
        else
            ysm(seg) = smooth(y(seg), win, 'moving');  % R2016b
        end
    end
end

