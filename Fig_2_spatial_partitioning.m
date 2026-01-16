
% Parameters
window_size = 10;  % Size of the sliding window in mm
overlap = 5;    % Overlap between windows in mm

% Initialize lists for x and y coordinates
all_med_coords = [];
all_trans_coords = [];
coefficient_pos = [];
coefficient_neg = [];

% minimum of ranges covered by electrodes from all subjects +/- overlap
x_min = -23;
x_max = 43;
y_min = -23;
y_max = 30;



% Generate sliding window centers with overlap
x_centers = x_min:overlap:x_max;
y_centers = y_min:overlap:y_max;

% Initialize matrix to store the number of Cluster 1 and total electrodes in each bin
pos_electrodes_in_bins = zeros(length(y_centers), length(x_centers));
neg_electrodes_in_bins = zeros(length(y_centers), length(x_centers));
total_electrodes_in_bins = zeros(length(y_centers), length(x_centers));


% Collect only subjects of interest
for i = 1:length(electrode_ID)
    n = subject_ID(i);
    e = electrode_ID(i);
    if subject(n).electrode(e).olf_coor>0 & ismember(n,[1:6])
    med_coor = subject(n).electrode(e).med_coor;
    trans_coor = subject(n).electrode(e).trans_coor;
    
    all_med_coords = [all_med_coords; med_coor];
    all_trans_coords = [all_trans_coords; trans_coor];
    coefficient_pos = [coefficient_pos; coeffs{1}(i)>0 ];
    coefficient_neg = [coefficient_neg; coeffs{1}(i)<0 ];
    else
    end
end


% Loop over sliding window centers to calculate the number of Cluster 1 and total electrodes
for x_idx = 1:length(x_centers)
    for y_idx = 1:length(y_centers)
        % Define the window range
        x_start = x_centers(x_idx) - window_size/2;
        x_end = x_centers(x_idx) + window_size/2;
        y_start = y_centers(y_idx) - window_size/2;
        y_end = y_centers(y_idx) + window_size/2;
        
        % Find electrodes within the current window
        in_window = all_med_coords >= x_start & all_med_coords < x_end & ...
                    all_trans_coords >= y_start & all_trans_coords < y_end;
        
        % Count total electrodes and positive/negative electrodes
        total_electrodes = sum(in_window);
        pos_electrodes = sum(in_window & coefficient_pos );
        neg_electrodes = sum(in_window & coefficient_neg );
        
        % Store the counts
        total_electrodes_in_bins(y_idx, x_idx) = total_electrodes;
        pos_electrodes_in_bins(y_idx, x_idx) = pos_electrodes;
        neg_electrodes_in_bins(y_idx, x_idx) = neg_electrodes;
    end
end

%% statistics and spatial partitioning for positive encoding cluster
% first calculate global null proportion
total_pos = sum(coefficient_pos);
total_electrodes_overall = length(coefficient_pos);
p0 = total_pos / total_electrodes_overall;

% minimum electrodes per bin to consider
Nmin = 1;

z_map = nan(size(total_electrodes_in_bins));
for x_idx = 1:length(x_centers)
    for y_idx = 1:length(y_centers)
        n = total_electrodes_in_bins(y_idx, x_idx);
        if n >= Nmin
            k = pos_electrodes_in_bins(y_idx, x_idx);
            mu = n * p0;
            sigma = sqrt(n * p0 * (1 - p0));
            if sigma > 0
                z_map(y_idx, x_idx) = (k - mu) / sigma;
            else
                z_map(y_idx, x_idx) = NaN;
            end
        end
    end
end

z_thresh = 1;  % ~p<0.025 one-sided

supra = (z_map > z_thresh);                 % supra-threshold bins
supra(isnan(z_map)) = false;                % ignore NaNs
supra(total_electrodes_in_bins < Nmin) = false;

% use 8-connected components
CC = bwconncomp(supra, 8);

cluster_masses_obs = zeros(CC.NumObjects,1);
for c = 1:CC.NumObjects
    idx = CC.PixelIdxList{c};
    cluster_masses_obs(c) = nansum(z_map(idx));
end

max_cluster_mass_obs = 0;
if ~isempty(cluster_masses_obs)
    max_cluster_mass_obs = max(cluster_masses_obs);
end
% If there is at least one cluster, find the largest (by cluster mass)
if ~isempty(cluster_masses_obs)
    % Index of largest cluster
    [~, maxClusterIdx] = max(cluster_masses_obs);
    cluster_idx = CC.PixelIdxList{maxClusterIdx};  % linear indices
    
    % Convert to row / col indices
    [row_idx, col_idx] = ind2sub(size(z_map), cluster_idx);

    % ---- Centroid ----
    centroid_x = mean(x_centers(col_idx));
    centroid_y = mean(y_centers(row_idx));

    % ---- Span in coordinate space ----
    x_min_cluster = min(x_centers(col_idx));
    x_max_cluster = max(x_centers(col_idx));
    y_min_cluster = min(y_centers(row_idx));
    y_max_cluster = max(y_centers(row_idx));

    fprintf('Largest positive cluster centroid: (x = %.2f, y = %.2f) mm\n', ...
        centroid_x, centroid_y);

    fprintf('Cluster span: x = [%.2f, %.2f] mm, y = [%.2f, %.2f] mm\n', ...
        x_min_cluster-overlap, x_max_cluster+overlap, y_min_cluster-overlap, y_max_cluster+overlap);
end
% After you’ve identified the winning cluster bins: row_idx, col_idx
% ---- Electrode-based "bulk" location of winning cluster ----
% Uses your existing arrays:
%   all_med_coords, all_trans_coords, coefficient_pos
% and the cluster bins:
%   row_idx, col_idx
% plus your window parameters:
%   window_size, x_centers, y_centers

% 1) Build a mask of bins belonging to the winning cluster
cluster_mask = false(size(z_map));
cluster_mask(sub2ind(size(z_map), row_idx, col_idx)) = true;

% 2) For each electrode, determine if it falls in ANY cluster bin (remember bins overlap)
in_cluster_any = false(size(all_med_coords));

for x_idx = 1:length(x_centers)
    for y_idx = 1:length(y_centers)
        if ~cluster_mask(y_idx, x_idx), continue; end

        x_start = x_centers(x_idx) - window_size/2;
        x_end   = x_centers(x_idx) + window_size/2;
        y_start = y_centers(y_idx) - window_size/2;
        y_end   = y_centers(y_idx) + window_size/2;

        in_window = all_med_coords  >= x_start & all_med_coords  < x_end & ...
                    all_trans_coords >= y_start & all_trans_coords < y_end;

        in_cluster_any = in_cluster_any | in_window;
    end
end

% 3) Decide whether "95% of electrodes" means ALL electrodes in the cluster footprint
%    or only POSITIVE electrodes within that footprint. Choose ONE:

use_pos_only = true;  % <-- set false if you want all electrodes in cluster bins

if use_pos_only
    keep = in_cluster_any & coefficient_pos;
else
    keep = in_cluster_any;
end

xc = all_med_coords(keep);
yc = all_trans_coords(keep);

fprintf('\nElectrodes in winning cluster footprint: %d\n', sum(in_cluster_any));
fprintf('Electrodes used for bulk summary: %d (use_pos_only=%d)\n', numel(xc), use_pos_only);

if numel(xc) < 5
    warning('Too few electrodes to define a stable 95%% region. Consider use_pos_only=false or lower threshold.');
else
    % 4) Robust center (median) and 95% core quantile box
    cx = median(xc, 'omitnan');
    cy = median(yc, 'omitnan');

    x95 = quantile(xc, [0.025 0.975]);
    y95 = quantile(yc, [0.025 0.975]);

    fprintf('Robust center (median): x=%.2f, y=%.2f mm\n', cx, cy);
    fprintf('95%% core (quantile box): x=[%.2f, %.2f] mm; y=[%.2f, %.2f] mm\n', ...
        x95(1), x95(2), y95(1), y95(2));
end




%% ---- Plot z-map with cluster outlines ----

z_plot = z_map;
alpha_mask = ~isnan(z_plot);   % 1 where data exists, 0 where background should show
%alpha_mask = z_plot > 1;
figure;
set(gcf, 'Color','w');   % figure background white

h = imagesc(x_centers, y_centers, z_plot);
set(h, 'AlphaData', alpha_mask);   % transparent NaNs

set(gca,'YDir','normal','Color','w');  % axes background white
colormap(jet);
clim([-2 2])
colorbar;

hold on;
xlabel('Medial/Lateral (mm)', 'FontSize', 14);
ylabel('Rostral/Caudal (mm)', 'FontSize', 14);
title('Binomial z-map with spatial clusters', 'FontSize', 14);
set(gca,'TickDir','out');
box off;


%% ---- Parameters for permutation test ----
Nmin    = 1;       % minimum electrodes per bin to consider
z_thresh = 1;   % threshold for supra-threshold bins
nPerm   = 2000;    % number of permutations

% Global null proportion of positive electrodes
total_pos = sum(coefficient_pos);
total_electrodes_overall = numel(coefficient_pos);
p0 = total_pos / total_electrodes_overall;

%% ---- Permutation loop ----
nElec = numel(coefficient_pos);
max_cluster_mass_perm = zeros(nPerm, 1);

for p = 1:nPerm
    
    % 1) Shuffle positive/negative labels across electrodes
    shuffled_pos = coefficient_pos(randperm(nElec));   % logical vector
    
    % 2) Recompute positive counts per bin under shuffled labels
    pos_electrodes_in_bins_perm = zeros(size(total_electrodes_in_bins));
    
    for x_idx = 1:length(x_centers)
        for y_idx = 1:length(y_centers)
            % same window definitions as original code
            x_start = x_centers(x_idx) - window_size/2;
            x_end   = x_centers(x_idx) + window_size/2;
            y_start = y_centers(y_idx) - window_size/2;
            y_end   = y_centers(y_idx) + window_size/2;
            
            % electrodes within current window (location only)
            in_window = all_med_coords >= x_start & all_med_coords < x_end & ...
                        all_trans_coords >= y_start & all_trans_coords < y_end;
            
            % number of positives in this bin under shuffled labels
            pos_electrodes_in_bins_perm(y_idx, x_idx) = sum(in_window & shuffled_pos);
            % NOTE: total_electrodes_in_bins is unchanged, so no need to recompute n
        end
    end
    
    % 3) Compute binomial z-map under permutation
    z_map_perm = nan(size(total_electrodes_in_bins));
    
    for x_idx = 1:length(x_centers)
        for y_idx = 1:length(y_centers)
            n = total_electrodes_in_bins(y_idx, x_idx);   % total count is fixed
            if n >= Nmin
                k = pos_electrodes_in_bins_perm(y_idx, x_idx);
                mu = n * p0;
                sigma = sqrt(n * p0 * (1 - p0));
                if sigma > 0
                    z_map_perm(y_idx, x_idx) = (k - mu) / sigma;
                end
            end
        end
    end
    
    % 4) Threshold and find clusters in permuted map
    supra_perm = (z_map_perm > z_thresh);
    supra_perm(isnan(z_map_perm)) = false;
    supra_perm(total_electrodes_in_bins < Nmin) = false;
    
    CC_perm = bwconncomp(supra_perm, 8);
    
    if CC_perm.NumObjects == 0
        max_cluster_mass_perm(p) = 0;
    else
        cluster_masses_perm = zeros(CC_perm.NumObjects,1);
        for c = 1:CC_perm.NumObjects
            idx = CC_perm.PixelIdxList{c};
            cluster_masses_perm(c) = nansum(z_map_perm(idx));
        end
        max_cluster_mass_perm(p) = max(cluster_masses_perm);
    end
end

% ---- Get p-value for largest observed cluster ----
p_global = mean(max_cluster_mass_perm >= max_cluster_mass_obs);
fprintf('Cluster-based permutation p (max cluster mass): %.4f\n', p_global);

%% statistics and spatial partitioning for negative encoding cluster
% first calculate global null proportion
total_neg = sum(coefficient_neg);
total_electrodes_overall = length(coefficient_neg);
p0 = total_neg / total_electrodes_overall;

% minimum electrodes per bin to consider
Nmin = 1;

z_map = nan(size(total_electrodes_in_bins));
for x_idx = 1:length(x_centers)
    for y_idx = 1:length(y_centers)
        n = total_electrodes_in_bins(y_idx, x_idx);
        if n >= Nmin
            k = neg_electrodes_in_bins(y_idx, x_idx);
            mu = n * p0;
            sigma = sqrt(n * p0 * (1 - p0));
            if sigma > 0
                z_map(y_idx, x_idx) = (k - mu) / sigma;
            else
                z_map(y_idx, x_idx) = NaN;
            end
        end
    end
end

z_thresh = 1;  % ~p<0.025 one-sided

supra = (z_map > z_thresh);                 % supra-threshold bins
supra(isnan(z_map)) = false;                % ignore NaNs
supra(total_electrodes_in_bins < Nmin) = false;

% use 8-connected components
CC = bwconncomp(supra, 8);

cluster_masses_obs = zeros(CC.NumObjects,1);
for c = 1:CC.NumObjects
    idx = CC.PixelIdxList{c};
    cluster_masses_obs(c) = nansum(z_map(idx));
end

max_cluster_mass_obs = 0;
if ~isempty(cluster_masses_obs)
    max_cluster_mass_obs = max(cluster_masses_obs);
end

% If there is at least one cluster, find the largest (by cluster mass)
if ~isempty(cluster_masses_obs)
    % Index of largest cluster
    [~, maxClusterIdx] = max(cluster_masses_obs);
    cluster_idx = CC.PixelIdxList{maxClusterIdx};  % linear indices
    
    % Convert to row / col indices
    [row_idx, col_idx] = ind2sub(size(z_map), cluster_idx);

    % ---- Centroid ----
    centroid_x = mean(x_centers(col_idx));
    centroid_y = mean(y_centers(row_idx));

    % ---- Span in coordinate space ----
    x_min_cluster = min(x_centers(col_idx));
    x_max_cluster = max(x_centers(col_idx));
    y_min_cluster = min(y_centers(row_idx));
    y_max_cluster = max(y_centers(row_idx));

    fprintf('Largest positive cluster centroid: (x = %.2f, y = %.2f) mm\n', ...
        centroid_x, centroid_y);

    fprintf('Cluster span: x = [%.2f, %.2f] mm, y = [%.2f, %.2f] mm\n', ...
        x_min_cluster-overlap, x_max_cluster+overlap, y_min_cluster-overlap, y_max_cluster+overlap);
end

%% ---- Plot z-map with cluster outlines ----

z_plot = z_map;
alpha_mask = ~isnan(z_plot);   % 1 where data exists, 0 where background should show
alpha_mask = z_plot>1;
figure;
set(gcf, 'Color','w');   % figure background white

h = imagesc(x_centers, y_centers, z_plot);
set(h, 'AlphaData', alpha_mask);   % transparent NaNs

set(gca,'YDir','normal','Color','w');  % axes background white
colormap(jet);
clim([-2 2])
colorbar;

hold on;
xlabel('Medial/Lateral (mm)', 'FontSize', 14);
ylabel('Rostral/Caudal (mm)', 'FontSize', 14);
title('Binomial z-map with spatial clusters', 'FontSize', 14);
set(gca,'TickDir','out');
box off;
%% ---- Parameters for permutation test ----
Nmin    = 1;       % minimum electrodes per bin to consider
z_thresh = 1;   % threshold for supra-threshold bins (one-sided ~p<0.025)
nPerm   = 2000;    % number of permutation

% Global null proportion of negitive electrodes
total_neg = sum(coefficient_neg);
total_electrodes_overall = numel(coefficient_neg);
p0 = total_neg / total_electrodes_overall;

%% ---- Permutation loop ----
nElec = numel(coefficient_neg);
max_cluster_mass_perm = zeros(nPerm, 1);

for p = 1:nPerm
    
    % 1) Shuffle negitive/negative labels across electrodes
    shuffled_neg = coefficient_neg(randperm(nElec));   % logical vector
    
    % 2) Recompute negitive counts per bin under shuffled labels
    neg_electrodes_in_bins_perm = zeros(size(total_electrodes_in_bins));
    
    for x_idx = 1:length(x_centers)
        for y_idx = 1:length(y_centers)
            % same window definitions as original code
            x_start = x_centers(x_idx) - window_size/2;
            x_end   = x_centers(x_idx) + window_size/2;
            y_start = y_centers(y_idx) - window_size/2;
            y_end   = y_centers(y_idx) + window_size/2;
            
            % electrodes within current window (location only)
            in_window = all_med_coords >= x_start & all_med_coords < x_end & ...
                        all_trans_coords >= y_start & all_trans_coords < y_end;
            
            % number of negitives in this bin under shuffled labels
            neg_electrodes_in_bins_perm(y_idx, x_idx) = sum(in_window & shuffled_neg);
            % NOTE: total_electrodes_in_bins is unchanged, so no need to recompute n
        end
    end
    
    % 3) Compute binomial z-map under permutation
    z_map_perm = nan(size(total_electrodes_in_bins));
    
    for x_idx = 1:length(x_centers)
        for y_idx = 1:length(y_centers)
            n = total_electrodes_in_bins(y_idx, x_idx);   % total count is fixed
            if n >= Nmin
                k = neg_electrodes_in_bins_perm(y_idx, x_idx);
                mu = n * p0;
                sigma = sqrt(n * p0 * (1 - p0));
                if sigma > 0
                    z_map_perm(y_idx, x_idx) = (k - mu) / sigma;
                end
            end
        end
    end
    
    % 4) Threshold and find clusters in permuted map
    supra_perm = (z_map_perm > z_thresh);
    supra_perm(isnan(z_map_perm)) = false;
    supra_perm(total_electrodes_in_bins < Nmin) = false;
    
    CC_perm = bwconncomp(supra_perm, 8);
    
    if CC_perm.NumObjects == 0
        max_cluster_mass_perm(p) = 0;
    else
        cluster_masses_perm = zeros(CC_perm.NumObjects,1);
        for c = 1:CC_perm.NumObjects
            idx = CC_perm.PixelIdxList{c};
            cluster_masses_perm(c) = nansum(z_map_perm(idx));
        end
        max_cluster_mass_perm(p) = max(cluster_masses_perm);
    end
end

% ---- Get p-value for largest observed cluster ----
p_global = mean(max_cluster_mass_perm >= max_cluster_mass_obs);
fprintf('Cluster-based permutation p (max cluster mass): %.4f\n', p_global);
