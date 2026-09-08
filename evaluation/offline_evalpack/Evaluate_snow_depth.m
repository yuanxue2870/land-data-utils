%% Unified Snow Depth Evaluation and Multi-Reference Scorecard Script
%% Author: Yuan Xue, Gemini, Github Copilot
%% Date: 05/08/2026

clc; clear;

%% Constants and fixed eval directories
LARGE_SNOW_LIMIT = 1e12;
land_eval_directory = '/scratch4/NCEPDEV/land/data/evaluation/';

% Output directory
output_dir = 'Images_snow_depth_eval';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

custom_colors = [1 1 1; jet(100)];
cmin = -50;
cmax = 50;
custom_colors_diff = b2r(cmin, cmax);
afs = 8; % fontsize in plots

%% User inputs
config_inputs; % Runs central config and loads all variables

% Map script-specific variations
output_step = eval_output_step;
exp_names   = eval_exp_names;

%% Derived variables
if ischar(exp_names)
    exp_names = strtrim(strsplit(exp_names, ','));
elseif isstring(exp_names)
    exp_names = cellstr(strtrim(exp_names));
end

num_experiments = numel(exp_names);
date_vector = start_date:output_step:finish_date;
total_days = numel(date_vector);

% Input Validations
if start_date > finish_date
    error('Start date must be before or equal to finish date.');
end
if ~strcmp(ocn_res, 'mx025')
    error('Only support mx025 for now');
end
if isnan(res) || ~ismember(res, [96, 384, 768, 1152])
    error('Invalid resolution. Must be one of: 96, 384, 768, 1152');
end
if isnan(output_step) || output_step <= 0
    error('Output step must be a positive number.');
end
if threshold_num_days < 1 || threshold_num_days > 60
    error('threshold_num_days must be between 1 and 60. Current value is %d.', threshold_num_days);
end
if threshold_num_days > total_days
    error('threshold_num_days (%d) cannot be greater than total number of days (%d).', threshold_num_days, total_days);
end
if ~exist(directory, 'dir')
    error('Experiment directory does not exist: %s', directory);
end
if ~exist(land_eval_directory, 'dir')
    error('Evaluation directory does not exist: %s', land_eval_directory);
end
if isempty(exp_names) || ~iscellstr(exp_names)
    error('Experiment names must be a non-empty cell array of strings.');
end

%% Determine vector length based on resolution
switch res
    case 1152
        vector_len = 2381853;
    case 96
        vector_len = 18320;
    case 384
        vector_len = 272699;
    case 768
        vector_len = 1067333;
    otherwise
        error('Unknown resolution: %d. Please check!', res);
end

%% Preallocate vector arrays
ex_sd_vector        = NaN(vector_len, total_days, num_experiments);
snodas_sd_vector    = NaN(vector_len, total_days);
era5_sd_vector      = NaN(vector_len, total_days);
era5_land_sd_vector = NaN(vector_len, total_days);

%% Load land fraction mask and coordinates
mask_dir = '/scratch3/NCEPDEV/stmp/Yuan.Xue/stats/';
filename = sprintf('C%d_hr3_%s_era5_input.mat', res, ocn_res);
mask_file = fullfile(mask_dir, filename);
if ~exist(mask_file, 'file')
    error('Mask file does not exist: %s', mask_file);
end
maskData = load(mask_file);
land_frac_mask = maskData.(sprintf('C%d_land_frac_FULL', res));
lat2D_corners  = maskData.(sprintf('C%d_lat2D_corners_FULL', res));
lon2D_corners  = maskData.(sprintf('C%d_lon2D_corners_FULL', res));
vtype_1D       = maskData.(sprintf('era5_C%d_vegetation_class', res));
lat1D          = maskData.(sprintf('C%d_lat1D', res));
lon1D          = maskData.(sprintf('C%d_lon1D', res));
num_tiles      = size(land_frac_mask, 3);

land_mask = ones(res, res, num_tiles);
ind_notOK = maskData.vegetation_class_FULL == 15 | maskData.vegetation_class_FULL <= 0;
land_mask(ind_notOK) = NaN;

land_mask_vector = ones(vector_len, 1);
ind_notOK_vec = vtype_1D == 15 | vtype_1D <= 0;
land_mask_vector(ind_notOK_vec) = NaN;

%% Load Bukovsky Mask Data (mask_data.nc)
bukovsky_nc_file = [mask_dir, 'mask_data.nc'];
if ~exist(bukovsky_nc_file, 'file')
    error('Bukovsky mask file mask_data.nc not found.');
end

buk_lon  = double(ncread(bukovsky_nc_file, 'lon'));
buk_lat  = double(ncread(bukovsky_nc_file, 'lat'));
buk_mask = ncread(bukovsky_nc_file, 'mask');
buk_attr = ncreadatt(bukovsky_nc_file, '/', 'region_list');

if ischar(buk_attr) || isstring(buk_attr)
    all_buk_regions = strtrim(strsplit(char(buk_attr), ','));
else
    all_buk_regions = buk_attr;
end


% Define land-centric Bukovsky basic regions
basic_regions = { ...
    'PacificNW',      'PacificSW',      'Southwest',      'Mezquital',     ...
    'NRockies',       'SRockies',       'GreatBasin',     'NPlains',       ...
    'CPlains',        'SPlains',        'Prairie',        'GreatLakes',    ...
    'Appalachia',     'DeepSouth',      'Southeast',      'MidAtlantic',   ...
    'EastBoreal',     'WestBoreal',     'EastTaiga',      'WestTaiga',     ...
    'CentralTundra',  'WestTundra',     'EastTundra'                       ...
    % Excluded: 'ColdNEPacific', 'WarmNEPacific', 'WarmNWAtlantic', 
    %           'NorthAtlantic', 'ColdNWAtlantic', 'Hudson', 'LabradorSea', 'Greenland'
};
x_regions = [{'NH', 'NA', 'CONUS'}, basic_regions];

% Map 1D Model Grid Points to Bukovsky 2D Grid
lon1D_norm = lon1D;
lon1D_norm(lon1D_norm < 0) = lon1D_norm(lon1D_norm < 0) + 360;
d_lat = buk_lat(2) - buk_lat(1);
d_lon = buk_lon(2) - buk_lon(1);

lat_idx_1D = round((lat1D - buk_lat(1)) / d_lat) + 1;
lon_idx_1D = round((lon1D_norm - buk_lon(1)) / d_lon) + 1;

valid_grid_pts = find(lat_idx_1D >= 1 & lat_idx_1D <= numel(buk_lat) & ...
                      lon_idx_1D >= 1 & lon_idx_1D <= numel(buk_lon) & ...
                      ~isnan(land_mask_vector));

lin_2d_idx = sub2ind([numel(buk_lon), numel(buk_lat)], ...
                     lon_idx_1D(valid_grid_pts), lat_idx_1D(valid_grid_pts));

fprintf('Starting run: res=%d, ocn_res=%s, total_days=%d, num_tiles=%d, experiments=%d\n', ...
        res, ocn_res, total_days, num_tiles, num_experiments);

%% Single-Pass Data Ingestion Loop
tic;
current_date = start_date;
counter = 0;

while current_date <= finish_date
    counter = counter + 1;
    fprintf('Processing day [%d/%d]: %s\n', counter, total_days, datestr(current_date, 'yyyy-mm-dd'));
    yyyys = datestr(current_date, 'yyyy');
    mms   = datestr(current_date, 'mm');
    dds   = datestr(current_date, 'dd');

    % Experiment files
    for mm = 1:num_experiments
        current_ex = exp_names{mm};
        file = sprintf('%s%s/mem000/restarts/vector/ufs_land_restart_back.%s-%s-%s_00-00-00.nc', ...
                      directory, current_ex, yyyys, mms, dds);
        if exist(file, 'file')
            try
                ex_sd = ncread(file, 'snow_depth');
                ex_sd(ex_sd > LARGE_SNOW_LIMIT | ex_sd < 0) = NaN;
                ex_sd_vector(:, counter, mm) = ex_sd;
            catch ME
                warning('Error reading experiment file %s: %s', file, ME.message);
            end
        else
            warning('Experiment file not found: %s', file);
        end
    end

    % SNODAS file
    snodas_file = sprintf('%s/SNODAS/C%d/SNODAS_snow.C%d.%s%s%s.nc', ...
                         land_eval_directory, res, res, yyyys, mms, dds);
    if exist(snodas_file, 'file')
        try
            snodas_sd = ncread(snodas_file, 'snow_depth');
            snodas_sd(snodas_sd < 0 | snodas_sd > LARGE_SNOW_LIMIT) = NaN;
            snodas_sd_vector(:, counter) = snodas_sd;
        catch ME
            warning('Error reading SNODAS file %s: %s', snodas_file, ME.message);
        end
    end

    % ERA5 file
    era5_file = sprintf('%s/ERA5_organize/snow/C%d.%s/C%d_era5_snow_%s%s%s00.nc', ...
                       land_eval_directory, res, ocn_res, res, yyyys, mms, dds);
    if exist(era5_file, 'file')
        try
            era5_swe = ncread(era5_file, 'SWE');
            era5_swe(era5_swe < 0 | era5_swe > LARGE_SNOW_LIMIT) = NaN;
            era5_snow_density = ncread(era5_file, 'SnowDensity');
            era5_sd_vector(:, counter) = era5_swe .* 1e6 ./ era5_snow_density;
        catch ME
            warning('Error reading ERA5 file %s: %s', era5_file, ME.message);
        end
    end

    % ERA5-Land file
    era5_land_file = sprintf('%s/ERA5-land_new/snow/C%d.%s/C%d_era5land_snow_%s%s%s00.nc', ...
                            land_eval_directory, res, ocn_res, res, yyyys, mms, dds);
    if exist(era5_land_file, 'file')
        try
            era5_land_swe = ncread(era5_land_file, 'SWE');
            era5_land_swe(era5_land_swe < 0 | era5_land_swe > LARGE_SNOW_LIMIT) = NaN;
            era5_land_snow_density = ncread(era5_land_file, 'SnowDensity');
            era5_land_sd_vector(:, counter) = era5_land_swe .* 1e6 ./ era5_land_snow_density;
        catch ME
            warning('Error reading ERA5-Land file %s: %s', era5_land_file, ME.message);
        end
    end

    current_date = current_date + output_step;
end
fprintf('Data reading completed in %.2f seconds.\n', toc);

%% Plot Time-Averaged Maps
fprintf('Plotting routines for temporal averages ...\n');
tic;

datasets = {'ERA5', 'ERA5land', 'SNODAS'};
mean_vecs = { ...
    mean(era5_sd_vector, 2, 'omitnan'), ...
    mean(era5_land_sd_vector, 2, 'omitnan'), ...
    mean(snodas_sd_vector, 2, 'omitnan') ...
};

for exp_idx = 1:num_experiments
    datasets{end+1} = sprintf('Experiment_%s', exp_names{exp_idx});
    mean_vecs{end+1} = mean(ex_sd_vector(:, :, exp_idx), 2, 'omitnan');
end

for ds_idx = 1:length(datasets)
    ds_name = datasets{ds_idx};
    current_vec = mean_vecs{ds_idx};
    tile_data = vec2tiles(current_vec, land_frac_mask, num_tiles, res);
    
    fprintf('  -> Plotting Time-averaged map: %s\n', ds_name);
    figure(1); clf;
    
    if strcmp(ds_name, 'SNODAS')
        usamap([30, 52], [-150, -60]);
    else
        worldmap([-90, 90], [-180, 180]);
    end
    
    current_matrix_full = NaN(res, res, num_tiles);
    for i = 1:num_tiles
        current_matrix = tile_data{i} .* squeeze(land_mask(:,:,i));
        current_matrix_full(:,:,i) = current_matrix;
        
        for j = 1:4
            geoshow(squeeze(lat2D_corners(j,:,:,i)), squeeze(lon2D_corners(j,:,:,i)), ...
                    current_matrix, 'DisplayType', 'texturemap');
        end
    end
   
    fct_finalize_map(custom_colors, [-0.001 500], 'Timeavg sd (mm)', afs, mask_dir);

    overall_mean = mean(current_matrix_full(:), 'omitnan');
    overall_min  = min(current_matrix_full(:), [], 'omitnan');
    overall_max  = max(current_matrix_full(:), [], 'omitnan');
    
    title(sprintf('%s time-averaged snow depth (mean=%.2f, min=%.2f, max=%.2f)', ...
                  strrep(ds_name, '_', ' '), overall_mean, overall_min, overall_max)); 
    
    set(gcf, 'PaperPositionMode', 'auto');
    print('-djpeg', '-r300', fullfile(output_dir, sprintf('Time_averaged_SD_%s.jpg', ds_name)));
    close(1);
end
fprintf('Temporal average map plotting completed in %.2f seconds.\n', toc);

%% Plot Difference Maps Between Experiments
fprintf('Calculating and plotting temporal average differences...\n');
tic;
diff_pairs = [];
for j = 2:num_experiments
    for i = 1:j-1
        diff_pairs = [diff_pairs; j, i];
    end
end

for p = 1:size(diff_pairs, 1)
    idx2 = diff_pairs(p, 1);
    idx1 = diff_pairs(p, 2);
    exp_name2 = exp_names{idx2};
    exp_name1 = exp_names{idx1};
    
    diff_vec = mean_vecs{3 + idx2} - mean_vecs{3 + idx1};
    diff_tiles = vec2tiles(diff_vec, land_frac_mask, num_tiles, res);
    
    figure(1); clf;
    worldmap([-90, 90], [-180, 180]);
    
    diff_matrix_full = NaN(res, res, num_tiles);
    for i = 1:num_tiles
        diff_tile = diff_tiles{i} .* squeeze(land_mask(:,:,i));
        diff_matrix_full(:,:,i) = diff_tile;
        diff_tile(isnan(diff_tile)) = 0;
        for j_corner = 1:4
            geoshow(squeeze(lat2D_corners(j_corner,:,:,i)), squeeze(lon2D_corners(j_corner,:,:,i)), ...
                    diff_tile, 'DisplayType', 'texturemap');
        end
    end
    
    fct_finalize_map(custom_colors_diff, [cmin cmax], 'Diff (mm)', afs, mask_dir);

    d_mean = mean(diff_matrix_full(:), 'omitnan');
    d_min  = min(diff_matrix_full(:), [], 'omitnan');
    d_max  = max(diff_matrix_full(:), [], 'omitnan');
    title(sprintf('Mean Diff: %s - %s\n(Mean=%.2f, Min=%.2f, Max=%.2f)', ...
                  exp_name2, exp_name1, d_mean, d_min, d_max));
    print('-djpeg', '-r300', fullfile(output_dir, sprintf('Diff_Avg_SD_%s_minus_%s.jpg', exp_name2, exp_name1)));
    close(1);
end
fprintf('Temporal average difference plotting completed in %.2f seconds.\n', toc);

%% Vectorized Metrics Calculation & Spatial Mapping
fprintf('Calculating spatial metrics maps...\n');
tic;

refs = {'ERA5', 'ERA5-Land', 'SNODAS'};
ref_vecs = {era5_sd_vector, era5_land_sd_vector, snodas_sd_vector};
metric_types = {'bias', 'mae', 'rmse', 'ubrmse', 'corr'};

for r = 1:numel(refs)
    for m = 1:num_experiments
        exp_name = exp_names{m};
        
        for mt = 1:numel(metric_types)
            cur_type = metric_types{mt};
            full_vector = compute_metric(ref_vecs{r}, ex_sd_vector(:,:,m), threshold_num_days, cur_type);
            
            metric_tiles = vec2tiles(full_vector, land_frac_mask, num_tiles, res);
            
            m_mean = mean(full_vector .* land_mask_vector, 'omitnan');
            m_min  = min(full_vector .* land_mask_vector, [], 'omitnan');
            m_max  = max(full_vector .* land_mask_vector, [], 'omitnan'); 
            
            fig_map = figure('Visible', 'off', 'Units', 'pixels', 'Position', [100 100 1000 600]); 
            clf;
            if strcmp(refs{r}, 'SNODAS'), usamap([30 52], [-150 -60]);
            else, worldmap([-90, 90], [-180, 180]); end
            
            cmap = [1 1 1; jet(100)]; 
            switch cur_type
                case 'bias'; c_lims = [-200 200];
                case 'corr'; c_lims = [0 1];
                otherwise;   c_lims = [0 150];
            end
            
            for i = 1:num_tiles
                tile_masked = metric_tiles{i} .* squeeze(land_mask(:,:,i));
                for j = 1:4
                    geoshow(squeeze(lat2D_corners(j,:,:,i)), squeeze(lon2D_corners(j,:,:,i)), ...
                            tile_masked, 'DisplayType', 'texturemap');
                end
            end
           
            fct_finalize_map(cmap, c_lims, [cur_type ' (mm)'], afs, mask_dir);

            title(sprintf('%s: %s vs %s\n(Mean: %.2f, Min: %.2f, Max: %.2f)', ...
                upper(cur_type), exp_name, refs{r}, m_mean, m_min, m_max));
            
            img_out = fullfile(output_dir, sprintf('Map_%s_%s_%s.jpg', cur_type, exp_name, refs{r}));
            print(fig_map, '-djpeg', '-r300', img_out);
            close(fig_map);
        end
    end
end
fprintf('Spatial metrics calculation and mapping completed in %.2f seconds.\n', toc);

%% RMSE Analysis by Vegetation Type (Snow-on Only: Ref > 0 & Exp > 0)
fprintf('Analyzing RMSE across vegetation types (snow-on only)...\n');
tic;
full_vegNames = {'Evergreen NL', 'Evergreen BL', 'Deciduous NL', 'Deciduous BL', ...
                 'Mixed Forest', 'Closed Shrub', 'Open Shrub', 'Woody Savannas', ...
                 'Savannas', 'Grass', 'Perm. Wetlands', 'Crop', 'Urban', ...
                 'Crop/Mosaic', 'Snow and Ice', 'Barren', 'Water', ...
                 'Wooded Tundra', 'Mixed Tundra', 'Bare Tundra'};

valid_veg_ids = vtype_1D(~isnan(land_mask_vector));
vegTypes = unique(valid_veg_ids(valid_veg_ids >= 1 & valid_veg_ids <= numel(full_vegNames)));
vegNames = full_vegNames(vegTypes);
nVeg = numel(vegTypes);

for r = 1:numel(refs)
    ref_name_clean = strrep(refs{r}, '-', '_');
    RMSEMean = NaN(nVeg, num_experiments);
    RMSEErr  = NaN(nVeg, num_experiments);
    ref_data_full = ref_vecs{r};
    
    for e = 1:num_experiments
        exp_data_full = ex_sd_vector(:, :, e);
        
        % 1. Mask dataset strictly for snow-on condition (ref > 0 & exp > 0)
        ref_snow = ref_data_full;
        exp_snow = exp_data_full;
        no_snow_mask = (ref_snow <= 0) | (exp_snow <= 0) | isnan(ref_snow) | isnan(exp_snow);
        ref_snow(no_snow_mask) = NaN;
        exp_snow(no_snow_mask) = NaN;
        
        % 2. Calculate cell temporal RMSE over snow-on days (min threshold = 1 day)
        cell_rmse = compute_metric(ref_snow, exp_snow, 1, 'rmse');
        
        % 3. Aggregate spatial statistics across vegetation types
        for v = 1:nVeg
            v_idx = (vtype_1D == vegTypes(v)) & ~isnan(land_mask_vector);
            
            % Extract valid cell RMSEs for this vegetation class
            v_cell_rmses = cell_rmse(v_idx);
            v_cell_rmses = v_cell_rmses(~isnan(v_cell_rmses));
            
            if ~isempty(v_cell_rmses)
                RMSEMean(v, e) = mean(v_cell_rmses);               % Spatial mean RMSE
                RMSEErr(v, e)  = std(v_cell_rmses) / sqrt(numel(v_cell_rmses)); % Spatial SEM across cells
            end
        end
    end

    figure('Units', 'pixels', 'Position', [100, 100, 1000, 600], 'Visible', 'off'); 
    b = bar(RMSEMean, 'grouped'); hold on;
    plot_colors = lines(num_experiments); 
    for e = 1:num_experiments
        b(e).FaceColor = plot_colors(e, :);
    end

    groupwidth = min(0.8, num_experiments / (num_experiments + 1.5));
    for e = 1:num_experiments
        x = (1:nVeg) - groupwidth/2 + (2*e-1)*groupwidth/(2*num_experiments);
        errorbar(x, RMSEMean(:,e), RMSEErr(:,e), 'k', 'linestyle', 'none', 'LineWidth', 1.0);
    end
    
    set(gca, 'XTick', 1:nVeg, 'XTickLabel', vegNames, 'XTickLabelRotation', 45);
    ylabel('RMSE (mm)');
    ylim([0 800]);
    title(sprintf('RMSE vs Vegetation Type (Ref > 0 & Exp > 0): %s', refs{r}));
    set(gca, 'Fontsize', 13);
    legend(exp_names, 'Location', 'NorthOutside', 'Orientation', 'horizontal');
    grid on;
    saveas(gcf, fullfile(output_dir, sprintf('Veg_Analysis_RMSE_%s.jpg', ref_name_clean)));
    close(gcf);
end
fprintf('Vegetation RMSE analysis completed in %.2f seconds.\n', toc);

%% Domain-Averaged Time Series Plots
fprintf('Generating domain-averaged time series plots...\n');
tic;
regions = {
    'North Dakota',      [45.9, 49.0], [-104.1, -96.5];
    'Montana',           [44.3, 49.0], [-116.0, -104.0];
    'Wyoming',           [41.0, 45.0], [-111.1, -104.0];
    'Colorado',          [37.0, 41.0], [-109.0, -102.0];
    'Alaska',            [54.0, 71.5], [-168.0, -130.0];
    'High mountain asia',[25.0, 45.0], [65.0, 105.0];
    'Turkey',            [36.0, 42.0], [26.0, 45.0];
    'Katzakastan',       [40.0, 56.0], [46.0, 87.0] 
};

for r_idx = 1:size(regions, 1)
    reg_name = regions{r_idx, 1};
    lat_lims = regions{r_idx, 2};
    lon_lims = regions{r_idx, 3};
    
    idx_in = find(lat1D >= lat_lims(1) & lat1D <= lat_lims(2) & ...
                  lon1D >= lon_lims(1) & lon1D <= lon_lims(2) & ...
                  ~isnan(land_mask_vector));
    
    if isempty(idx_in), continue; end
    
    figure('Visible', 'off', 'Units', 'pixels', 'Position', [100 100 1000 600]); hold on;
    plot(date_vector, mean(era5_sd_vector(idx_in, :), 1, 'omitnan'), 'k-', 'LineWidth', 2, 'DisplayName', 'ERA5');
    plot(date_vector, mean(era5_land_sd_vector(idx_in, :), 1, 'omitnan'), 'k--', 'LineWidth', 1.5, 'DisplayName', 'ERA5-Land');
    plot(date_vector, mean(snodas_sd_vector(idx_in, :), 1, 'omitnan'), 'r-', 'LineWidth', 2, 'DisplayName', 'SNODAS');
    
    exp_colors = lines(num_experiments);
    for m = 1:num_experiments
        plot(date_vector, mean(ex_sd_vector(idx_in, :, m), 1, 'omitnan'), 'Color', exp_colors(m,:), 'LineWidth', 1.5, 'DisplayName', exp_names{m});
    end
    
    datetick('x', 'mm/dd', 'keepticks');
    title(['Domain avg. time Series: ', reg_name]);
    ylabel('Snow Depth (mm)');
    legend('Location', 'northeastoutside');
    grid on;
    set(gca, 'fontsize', 13);
    saveas(gcf, fullfile(output_dir, sprintf('Timeseries_SD_%s.jpg', strrep(reg_name, ' ', '_'))));
    close(gcf);
end
fprintf('Domain-averaged time series plotting completed in %.2f seconds.\n', toc);

%% Multi-Reference Bukovsky Region Scorecards
fprintf('Generating Multi-Reference Bukovsky Regional Scorecards...\n');
tic;
color_dark_green  = [0.00, 0.50, 0.00];
color_light_green = [0.60, 0.95, 0.60];
color_white       = [1.00, 1.00, 1.00];
color_light_red   = [1.00, 0.75, 0.78];
color_dark_red    = [0.85, 0.00, 0.00];
color_gray        = [0.82, 0.82, 0.82];

scorecard_ref_names = {'ERA5-Land', 'SNODAS'};
scorecard_ref_vecs  = {era5_land_sd_vector, snodas_sd_vector};
metrics_list        = {'rmse', 'bias', 'mae', 'ubrmse', 'corr'};

baseline_exp = exp_names{1};
comp_exp     = exp_names{end};

for m = 1:numel(metrics_list)
    metric_type = metrics_list{m};
    metric_title = ternary(strcmp(metric_type, 'bias'), '|bias|', upper(metric_type));
    lower_is_better = ~strcmp(metric_type, 'corr');
    
    fig_w = 1250; fig_h = 750;
    fig_grid = figure('Units', 'pixels', 'Position', [100, 100, fig_w, fig_h], 'Visible', 'off', 'Color', 'w');
    clf; hold on; axis off; xlim([0, fig_w]); ylim([0, fig_h]);

    x_start = 140; y_start_top = 540; cell_w = 30; cell_h = 35;
    text(fig_w/2, fig_h - 40, sprintf('Scorecard for Snow Depth (%s)', metric_title), 'HorizontalAlignment', 'center', 'FontSize', 16, 'FontWeight', 'bold');
    text(fig_w/2, fig_h - 70, '(Lead Day = 00)', 'HorizontalAlignment', 'center', 'FontSize', 13);

    snodas_incomplete_regions = false(1, numel(x_regions));

    for k = 1:numel(scorecard_ref_names)
        cur_ref_name = scorecard_ref_names{k};
        cur_ref_vec  = scorecard_ref_vecs{k};
        cy = y_start_top - (k - 1) * cell_h;
        text(x_start - 10, cy + cell_h/2, cur_ref_name, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', 'FontSize', 10, 'FontWeight', 'bold');
        
        for r = 1:numel(x_regions)
            reg_name = x_regions{r};
            cx = x_start + (r - 1) * cell_w;
            
            if strcmp(cur_ref_name, 'SNODAS') && ismember(reg_name, {'NH', 'NA'})
                vec_idx = [];
            else
                switch reg_name
                    case 'NH'
                        vec_idx = valid_grid_pts(lat1D(valid_grid_pts) >= 0);
                    case 'NA'
                        vec_idx = valid_grid_pts(lat1D(valid_grid_pts) >= 15 & lat1D(valid_grid_pts) <= 85 & lon1D_norm(valid_grid_pts) >= 190 & lon1D_norm(valid_grid_pts) <= 310);
                    case 'CONUS'
                        vec_idx = valid_grid_pts(lat1D(valid_grid_pts) >= 24.5 & lat1D(valid_grid_pts) <= 49.5 & lon1D_norm(valid_grid_pts) >= 235 & lon1D_norm(valid_grid_pts) <= 293);
                        has_snodas_data = any(~isnan(snodas_sd_vector(vec_idx, :)), 2);
                        vec_idx = vec_idx(has_snodas_data);
                    otherwise
                        layer_k = find(strcmpi(all_buk_regions, reg_name), 1);
                        if ~isempty(layer_k)
                            mask_layer_2d = buk_mask(:, :, layer_k);
                            vec_idx = valid_grid_pts(mask_layer_2d(lin_2d_idx) == 1);
                        else
                            vec_idx = [];
                        end
                end
            end
            
            if isempty(vec_idx)
                bg_col = color_gray; sig_level = 0; is_better = false;
                if strcmp(cur_ref_name, 'SNODAS'), snodas_incomplete_regions(r) = true; end
            else
                base_val = compute_metric(cur_ref_vec(vec_idx, :), ex_sd_vector(vec_idx, :, 1), threshold_num_days, metric_type);
                comp_val = compute_metric(cur_ref_vec(vec_idx, :), ex_sd_vector(vec_idx, :, end), threshold_num_days, metric_type);
                
                valid  = ~isnan(base_val) & ~isnan(comp_val);
                data_b = base_val(valid); data_c = comp_val(valid);
                
                is_100pct_covered = ternary(strcmp(cur_ref_name, 'SNODAS'), all(any(~isnan(cur_ref_vec(vec_idx, :)), 2)), true);
                if ~is_100pct_covered, snodas_incomplete_regions(r) = true; end

                if isempty(data_c) || length(data_c) < 2 || ~is_100pct_covered
                    bg_col = color_gray; sig_level = 0; is_better = false;
                else
                    delta = mean(data_c) - mean(data_b);
                    is_better = ternary(strcmp(metric_type, 'bias'), abs(mean(data_c)) < abs(mean(data_b)), ternary(lower_is_better, delta < 0, delta > 0));
                    
                    [~, p_val] = ttest(data_c, data_b, 'Alpha', 0.05);
                    sig_level = ternary(isnan(p_val), 0, ternary(p_val < 0.001, 3, ternary(p_val < 0.01, 2, ternary(p_val < 0.05, 1, 0))));
                    bg_col = ternary(sig_level == 1, ternary(is_better, color_light_green, color_light_red), ternary(sig_level >= 2, color_white, color_gray));
                end
            end
            
            rectangle('Position', [cx, cy, cell_w, cell_h], 'FaceColor', bg_col, 'EdgeColor', 'k', 'LineWidth', 1.1);
            mid_x = cx + cell_w / 2; mid_y = cy + cell_h / 2;
            
            if sig_level >= 2
                sym_color = ternary(is_better, color_dark_green, color_dark_red);
                if is_better
                    % Upright triangle (Improvement)
                    if sig_level == 3
                        patch([mid_x - 7, mid_x + 7, mid_x], [mid_y - 7, mid_y - 7, mid_y + 8], sym_color, 'EdgeColor', 'none');
                    elseif sig_level == 2
                        patch([mid_x - 4.5, mid_x + 4.5, mid_x], [mid_y - 4.5, mid_y - 4.5, mid_y + 5.5], sym_color, 'EdgeColor', 'none');
                    end
                else
                    % Upside-down triangle (Degradation)
                    if sig_level == 3
                        patch([mid_x - 7, mid_x + 7, mid_x], [mid_y + 7, mid_y + 7, mid_y - 8], sym_color, 'EdgeColor', 'none');
                    elseif sig_level == 2
                        patch([mid_x - 4.5, mid_x + 4.5, mid_x], [mid_y + 4.5, mid_y + 4.5, mid_y - 5.5], sym_color, 'EdgeColor', 'none');
                    end
                end
            end

            if k == numel(scorecard_ref_names)
                display_label = reg_name;
                if snodas_incomplete_regions(r), display_label = [reg_name, '*']; end
                text(cx + cell_w/2, cy - 12, display_label, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'Rotation', 60, 'FontSize', 9, 'FontWeight', 'bold');
            end
        end
    end

    rectangle('Position', [x_start, y_start_top - (numel(scorecard_ref_names) - 1) * cell_h, numel(x_regions) * cell_w, numel(scorecard_ref_names) * cell_h], 'EdgeColor', 'k', 'LineWidth', 1.5);
    line([x_start + 3 * cell_w, x_start + 3 * cell_w], [y_start_top - (numel(scorecard_ref_names) - 1) * cell_h, y_start_top + cell_h], 'Color', 'k', 'LineWidth', 2);

    %% DRAW VISUAL LEGEND
    leg_y_top = 210;
    leg_x_left = 80;

    legend_items = {
        color_white,       'tri_lg_g', sprintf('%s is better than %s at the 99.9%% significance level', comp_exp, baseline_exp);
        color_white,       'tri_sm_g', sprintf('%s is better than %s at the 99%% significance level', comp_exp, baseline_exp);
        color_light_green, 'bg_only',  sprintf('%s is better than %s at the 95%% significance level', comp_exp, baseline_exp);
        color_white,       'tri_lg_r', sprintf('%s is worse than %s at the 99.9%% significance level', comp_exp, baseline_exp);
        color_white,       'tri_sm_r', sprintf('%s is worse than %s at the 99%% significance level', comp_exp, baseline_exp);
        color_light_red,   'bg_only',  sprintf('%s is worse than %s at the 95%% significance level', comp_exp, baseline_exp);
        color_gray,        'bg_only',  sprintf('No statistically significant difference / Incomplete regional coverage')
    };

    for k = 1:size(legend_items, 1)
        iy = leg_y_top - (k - 1) * 22;
        bg_c   = legend_items{k, 1};
        style  = legend_items{k, 2};
        txt_label = legend_items{k, 3};
        
        rectangle('Position', [leg_x_left, iy, 16, 16], 'FaceColor', bg_c, 'EdgeColor', 'k', 'LineWidth', 1);
        
        lx = leg_x_left + 8;
        ly = iy + 8;
        
        if strcmp(style, 'tri_lg_g')
            patch([lx - 6, lx + 6, lx], [ly - 6, ly - 6, ly + 7], color_dark_green, 'EdgeColor', 'none');
        elseif strcmp(style, 'tri_sm_g')
            patch([lx - 4, lx + 4, lx], [ly - 4, ly - 4, ly + 5], color_dark_green, 'EdgeColor', 'none');
        elseif strcmp(style, 'tri_lg_r')
            patch([lx - 6, lx + 6, lx], [ly + 6, ly + 6, ly - 7], color_dark_red, 'EdgeColor', 'none');
        elseif strcmp(style, 'tri_sm_r')
            patch([lx - 4, lx + 4, lx], [ly + 4, ly + 4, ly - 5], color_dark_red, 'EdgeColor', 'none');
        end
        
        text(leg_x_left + 24, ly, txt_label, 'VerticalAlignment', 'middle', 'FontSize', 9, 'Interpreter', 'none');
    end

    % Legend Footnote for Asterisk
    note_y = leg_y_top - size(legend_items, 1) * 22 - 6;
    text(leg_x_left, note_y, '* Indicates region where SNODAS does not have 100% spatial grid coverage (SNODAS evaluation omitted/skipped).', ...
        'VerticalAlignment', 'top', 'FontSize', 8.5, 'FontAngle', 'italic', 'Interpreter', 'none');

    grid_filename = fullfile(output_dir, sprintf('Scorecard_LeadDay0_MultiRef_%s.jpg', upper(metric_type)));
    print(fig_grid, '-djpeg', '-r300', grid_filename); 
    close(fig_grid);
end
fprintf('Multi-reference regional scorecards completed in %.2f seconds.\n', toc);

%% Consolidated HTML Dashboard Generation
fprintf('Generating Master Dashboard HTML...\n');
tic;
dashboard_file = fullfile(output_dir, 'snow_depth_evaluation_dashboard.html');
fid = fopen(dashboard_file, 'w');

fprintf(fid, '<html><head><title>Snow Depth Evaluation Dashboard</title>\n');
fprintf(fid, '<style>\n');
fprintf(fid, '  body { font-family: Arial, sans-serif; margin: 20px; background-color: #f4f4f9; }\n');
fprintf(fid, '  h1, h2 { color: #333; text-align: center; }\n');
fprintf(fid, '  .section { background: #fff; padding: 20px; margin-bottom: 30px; border-radius: 8px; box-shadow: 0 2px 5px rgba(0,0,0,0.1); }\n');
fprintf(fid, '  .grid { display: flex; flex-wrap: wrap; justify-content: center; gap: 10px; }\n');
fprintf(fid, '  .card { border: 1px solid #ddd; padding: 10px; text-align: center; background: #fff; width: 400px; }\n');
fprintf(fid, '  img { width: 100%%; height: auto; border-radius: 4px; }\n');
fprintf(fid, '  .card-title { font-weight: bold; margin-bottom: 8px; font-size: 14px; color: #555; }\n');
fprintf(fid, '</style></head><body>\n');

fprintf(fid, '<h1>Snow Depth Evaluation: Offline land DA workflow at C%d</h1>\n', res);
fprintf(fid, '<p style="text-align:center">Date range: %s - %s</p>\n', datestr(start_date, 'yyyy-mm-dd'), datestr(finish_date, 'yyyy-mm-dd'));
fprintf(fid, '<p style="text-align:center">Generated on: %s</p>\n', datestr(now));

% 1. Temporal Averages
fprintf(fid, '<div class="section"><h2>Temporal Averages</h2><div class="grid">\n');
for ds_idx = 1:length(datasets)
    img_name = sprintf('Time_averaged_SD_%s.jpg', datasets{ds_idx});
    actual_img = get_existing_filename(img_name, output_dir);
    if ~isempty(actual_img)
        fprintf(fid, '  <div class="card"><div class="card-title">%s</div><img src="%s"></div>\n', datasets{ds_idx}, actual_img);
    end
end
fprintf(fid, '</div></div>\n');

% 2. Difference Maps
fprintf(fid, '<div class="section"><h2>Temporal Average Differences</h2><div class="grid">\n');
for p = 1:size(diff_pairs, 1)
    img_name = sprintf('Diff_Avg_SD_%s_minus_%s.jpg', exp_names{diff_pairs(p,1)}, exp_names{diff_pairs(p,2)});
    actual_img = get_existing_filename(img_name, output_dir);
    if ~isempty(actual_img)
        fprintf(fid, '<div class="card"><div class="card-title">%s minus %s</div><img src="%s"></div>\n', exp_names{diff_pairs(p,1)}, exp_names{diff_pairs(p,2)}, actual_img);
    end
end
fprintf(fid, '</div></div>\n');

% 3. Multi-Reference Regional Scorecards
fprintf(fid, '<div class="section"><h2>Multi-Reference Regional Scorecards</h2><div class="grid">\n');
for m = 1:numel(metrics_list)
    img_name = sprintf('Scorecard_LeadDay0_MultiRef_%s.jpg', upper(metrics_list{m}));
    actual_img = get_existing_filename(img_name, output_dir);
    if ~isempty(actual_img)
        fprintf(fid, '  <div class="card" style="width:600px;"><div class="card-title">%s Scorecard</div><img src="%s"></div>\n', upper(metrics_list{m}), actual_img);
    end
end
fprintf(fid, '</div></div>\n');

% 4. Spatial Metric Maps
for r = 1:numel(refs)
    fprintf(fid, '<div class="section"><h2>Metrics vs %s</h2><div class="grid">\n', refs{r});
    for m = 1:num_experiments
        for mt = 1:numel(metric_types)
            img_name = sprintf('Map_%s_%s_%s.jpg', metric_types{mt}, exp_names{m}, refs{r});
            actual_img = get_existing_filename(img_name, output_dir);
            if ~isempty(actual_img)
                fprintf(fid, '  <div class="card"><div class="card-title">%s: %s vs %s</div><img src="%s"></div>\n', upper(metric_types{mt}), exp_names{m}, refs{r}, actual_img);
            end
        end
    end
    fprintf(fid, '</div></div>\n');
end

% 5. Vegetation Analysis
fprintf(fid, '<div class="section"><h2>RMSE by Vegetation Type</h2><div class="grid">\n');
for r = 1:numel(refs)
    img_name = sprintf('Veg_Analysis_RMSE_%s.jpg', strrep(refs{r}, '-', '_'));
    actual_img = get_existing_filename(img_name, output_dir);
    if ~isempty(actual_img)
        fprintf(fid, '  <div class="card" style="width:500px;"><div class="card-title">Reference: %s</div><img src="%s"></div>\n', refs{r}, actual_img);
    end
end
fprintf(fid, '</div></div>\n');

% 6. Time Series
fprintf(fid, '<div class="section"><h2>Regional Time Series</h2><div class="grid">\n');
for r_idx = 1:size(regions, 1)
    img_name = sprintf('Timeseries_SD_%s.jpg', strrep(regions{r_idx, 1}, ' ', '_'));
    actual_img = get_existing_filename(img_name, output_dir);
    if ~isempty(actual_img)
        fprintf(fid, '  <div class="card" style="width:450px;"><div class="card-title">%s</div><img src="%s"></div>\n', regions{r_idx, 1}, actual_img);
    end
end
fprintf(fid, '</div></div>\n');

fprintf(fid, '</body></html>\n');
fclose(fid);
fprintf('Master dashboard generated in %.2f seconds: %s\nAll tasks finished successfully!\n', toc, dashboard_file);


%% =========================================================================
%% HELPER FUNCTIONS
%% =========================================================================

function cmap = b2r(cmin, cmax)
    % B2R Generates a symmetric/asymmetric Blue-to-Red colormap centered at zero.
    if cmin >= cmax
        error('b2r:InvalidLimits', 'cmin must be less than cmax.');
    end
    
    num_colors = 256;
    
    if cmin >= 0
        r = linspace(1, 1, num_colors)';
        g = linspace(1, 0, num_colors)';
        b = linspace(1, 0, num_colors)';
        cmap = [r g b];
        return;
    elseif cmax <= 0
        r = linspace(0, 1, num_colors)';
        g = linspace(0, 1, num_colors)';
        b = linspace(1, 1, num_colors)';
        cmap = [r g b];
        return;
    end
    
    ratio = abs(cmin) / (cmax - cmin);
    num_blue = round(num_colors * ratio);
    num_red  = num_colors - num_blue;
    
    % Blue to White
    r_blue = linspace(0, 1, num_blue)';
    g_blue = linspace(0, 1, num_blue)';
    b_blue = linspace(1, 1, num_blue)';
    
    % White to Red
    r_red = linspace(1, 1, num_red)';
    g_red = linspace(1, 0, num_red)';
    b_red = linspace(1, 0, num_red)';
    
    cmap = [r_blue, g_blue, b_blue; r_red, g_red, b_red];
end

function actual_filename = get_existing_filename(target_filename, output_dir)
    full_path = fullfile(output_dir, target_filename);
    if exist(full_path, 'file')
        actual_filename = target_filename;
        return;
    end
    [~, name, ext] = fileparts(target_filename);
    search_pattern = [name, ext];
    d = dir(fullfile(output_dir, ['*', ext]));
    for k = 1:numel(d)
        if strcmpi(d(k).name, search_pattern)
            actual_filename = d(k).name;
            return;
        end
    end
    actual_filename = ''; % Not found
end

function tiles = vec2tiles(vec, land_frac_mask, num_tiles, res)
    tiles = cell(1, num_tiles);
    iloc = 0;
    for i = 1:num_tiles
        valid_mask_2D = land_frac_mask(:, :, i) > 0;
        valid_mask_1D = valid_mask_2D(:);
        num_valid = sum(valid_mask_1D);
        
        tile_data = NaN(res * res, 1);
        if num_valid > 0
            tile_data(valid_mask_1D) = vec(iloc + 1 : iloc + num_valid);
            iloc = iloc + num_valid;
        end
        tiles{i} = reshape(tile_data, [res, res]);
    end
end

function result = compute_metric(obs, exp, thresh, type)
    valid = ~isnan(obs) & ~isnan(exp);
    valid_cnt = sum(valid, 2);
    
    obs_clean = obs; obs_clean(~valid) = 0;
    exp_clean = exp; exp_clean(~valid) = 0;
    diff_v = exp_clean - obs_clean;
    
    switch lower(type)
        case 'bias'
            result = sum(diff_v, 2) ./ valid_cnt;
        case 'mae'
            result = sum(abs(diff_v), 2) ./ valid_cnt;
        case 'rmse'
            result = sqrt(sum(diff_v.^2, 2) ./ valid_cnt);
        case 'ubrmse'
            bias = sum(diff_v, 2) ./ valid_cnt;
            rmse2 = sum(diff_v.^2, 2) ./ valid_cnt;
            result = sqrt(max(0, rmse2 - bias.^2));
        case 'corr'
            mean_obs = sum(obs_clean, 2) ./ valid_cnt;
            mean_exp = sum(exp_clean, 2) ./ valid_cnt;
            
            dev_obs = (obs - mean_obs); dev_obs(~valid) = 0;
            dev_exp = (exp - mean_exp); dev_exp(~valid) = 0;
            
            cov_val = sum(dev_obs .* dev_exp, 2);
            var_obs = sum(dev_obs.^2, 2);
            var_exp = sum(dev_exp.^2, 2);
            
            denom = sqrt(var_obs .* var_exp);
            result = cov_val ./ denom;
            result(denom == 0) = NaN;
    end
    
    result(valid_cnt < thresh) = NaN;
end

function val = ternary(cond, true_val, false_val)
    if cond
        val = true_val;
    else
        val = false_val;
    end
end

function fct_finalize_map(custom_colors, colorbar_limits, colorbar_string, afs, mask_dir)
    persistent is_cached cached_land_lat cached_land_lon ...
               cached_intl_lat cached_intl_lon ...
               cached_usa_lat  cached_usa_lon

    if isempty(is_cached)
        if nargin < 5, mask_dir = ''; end
        [cached_land_lat, cached_land_lon] = load_shp_coords(mask_dir, 'landareas');
        [cached_intl_lat, cached_intl_lon] = load_shp_coords(mask_dir, 'intl_boundaries');
        [cached_usa_lat,  cached_usa_lon]  = load_shp_coords(mask_dir, 'usastatelo');
        is_cached = true;
    end

    set(gcf, 'Colormap', custom_colors);
    hl = colorbar('Location', 'EastOutside', 'FontSize', afs);
    clim(colorbar_limits);
    set(get(hl, 'ylabel'), 'String', colorbar_string, 'FontSize', afs);

    setm(gca, 'LabelFormat', 'signed', 'FontSize', afs);
    if ~isempty(cached_land_lat), geoshow(cached_land_lat, cached_land_lon, 'Color', 'black'); end
    if ~isempty(cached_intl_lat), geoshow(cached_intl_lat, cached_intl_lon, 'Color', 'black', 'LineWidth', 0.5); end
    if ~isempty(cached_usa_lat),  geoshow(cached_usa_lat,  cached_usa_lon,  'Color', 'black', 'LineWidth', 0.5); end
end

function [lat, lon] = load_shp_coords(mask_dir, filename)
    % Look in mask_dir first, fallback to MATLAB system path
    shp_path = fullfile(mask_dir, [filename, '.shp']);
    if ~exist(shp_path, 'file')
        shp_path = filename;
    end

    try
        S = shaperead(shp_path, 'UseGeoCoords', true);
        lat = [S.Lat];
        lon = [S.Lon];
    catch ME
        warning('Failed to load shapefile %s: %s', shp_path, ME.message);
        lat = []; lon = [];
    end
end
