%% The script is used to analyze diag files
%% Author: Yuan Xue, Gemini, Github Copilot
%% Date: 05/11/2026
%% Note: Gemini and Github Copilot were used to assist with developing this code. The code has been reviewed, edited, and validated by NWS staff.

clc; clear;

%% Constants and fixed eval directories
LARGE_SNOW_LIMIT = 1e12;
afs = 10;
cmin = 0;
cmax = 10;
custom_colors = jet(100);

% Output directory
output_dir = 'Images_diag_analysis';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

%% User inputs
config_inputs; % Runs the central config and loads all variables

% Map script-specific variations
output_step = diag_output_step;
exp_names   = eval_exp_names;

%% Derived variables
if ischar(exp_names)
    exp_names = strtrim(strsplit(exp_names, ','));
elseif isstring(exp_names)
    exp_names = cellstr(strtrim(exp_names));
end

num_experiments = numel(exp_names);
total_steps = numel(start_date:output_step:finish_date);
date_vector = start_date:output_step:finish_date;
date_range_str = sprintf('%s-%s', datestr(start_date, 'yyyymmdd'), datestr(finish_date, 'yyyymmdd'));

% Input Validations
if start_date > finish_date
    error('Start date must be before or equal to finish date.');
end
if ~strcmp(ocn_res, 'mx025')
    error('Only support mx025 for now');
end
if isnan(res) || ~ismember(res, [384, 768, 1152])
    error('Invalid resolution. Only support high-res. Must be one of: 384, 768, 1152');
end

%% Load land fraction mask
mask_dir = '/scratch3/NCEPDEV/stmp/Yuan.Xue/stats/';
mask_file = sprintf('%sC%d_hr3_%s_era5_input.mat', mask_dir, res, ocn_res);
if ~exist(mask_file, 'file')
    error('Mask file does not exist: %s', mask_file);
end
maskData = load(mask_file);
lat1D = maskData.(sprintf('C%d_lat1D', res));
lon1D = maskData.(sprintf('C%d_lon1D', res));

fprintf('Starting run: res=%d, total_steps=%d, experiments=%d\n', ...
        res, total_steps, num_experiments);

%% Process each date
tic; 
current_date = start_date;
counter = 0;

% Initialize data cubes: [GridPoints x Sources x Exps x Time]
oman_beforeQC_vector_FULL = NaN(length(lat1D), length(diag_source), num_experiments, total_steps);
ombg_beforeQC_vector_FULL = NaN(length(lat1D), length(diag_source), num_experiments, total_steps);
oman_afterQC_vector_FULL  = NaN(length(lat1D), length(diag_source), num_experiments, total_steps);
ombg_afterQC_vector_FULL  = NaN(length(lat1D), length(diag_source), num_experiments, total_steps);

% Initialize counts: [1 x Sources x Exps x Time]
counts_beforeQC = NaN(1, length(diag_source), num_experiments, total_steps);
counts_afterQC  = NaN(1, length(diag_source), num_experiments, total_steps);

% Pre-build KD-tree for the target grid
Mdl = KDTreeSearcher([lat1D, lon1D]);

while current_date <= finish_date
    counter = counter + 1;
    fprintf('Processing date: %s\n', datestr(current_date));
    
    yyyys = datestr(current_date, 'yyyy');
    mms   = datestr(current_date, 'mm');
    dds   = datestr(current_date, 'dd');
    hhs   = datestr(current_date, 'HH');

    for mm = 1:num_experiments
        current_ex = exp_names{mm};
        for nn = 1:length(diag_source)
            current_diag = diag_source{nn};
            
            file = sprintf('%s%s/DA/jedi_anl/diag_%s_%s%s%s%s.nc', ...
                          directory, current_ex, current_diag, yyyys, mms, dds, hhs);
            
            if exist(file, 'file')
                % 1. Read Data
                oman = ncread(file, 'oman/totalSnowDepth');
                ombg = ncread(file, 'ombg/totalSnowDepth');
                lat  = ncread(file, 'MetaData/latitude');
                lon  = ncread(file, 'MetaData/longitude');
                lon  = mod((lon+180),360)-180;
                
                % Outlier Filter
                oman(abs(oman) > LARGE_SNOW_LIMIT) = NaN;
                ombg(abs(ombg) > LARGE_SNOW_LIMIT) = NaN;
                
                % 2. Map Before QC (Double Averaging)
                [uni_locs_all, ~, ic_all] = unique([lat, lon], 'rows');
                oman_uni_all = accumarray(ic_all, oman, [], @(x) mean(x, 'omitnan'));
                ombg_uni_all = accumarray(ic_all, ombg, [], @(x) mean(x, 'omitnan'));
                counts_beforeQC(1,nn,mm,counter) = length(uni_locs_all);

                % Level 2: Average by Grid Index
                ind_min_all = knnsearch(Mdl, uni_locs_all);
                [grid_idx_all, ~, ic_grid_all] = unique(ind_min_all);
                oman_grid_all = accumarray(ic_grid_all, oman_uni_all, [], @(x) mean(x, 'omitnan'));
                ombg_grid_all = accumarray(ic_grid_all, ombg_uni_all, [], @(x) mean(x, 'omitnan'));

                oman_beforeQC_vector_FULL(grid_idx_all, nn, mm, counter) = oman_grid_all;
                ombg_beforeQC_vector_FULL(grid_idx_all, nn, mm, counter) = ombg_grid_all;

                % 3. Apply QC 
                qc0 = ncread(file, 'EffectiveQC0/totalSnowDepth');
                qc1 = ncread(file, 'EffectiveQC1/totalSnowDepth');
                indOK = find(qc0==0 & qc1==0);

                if ~isempty(indOK)
                    lat_ok = lat(indOK); lon_ok = lon(indOK);
                    oman_ok = oman(indOK); ombg_ok = ombg(indOK);

                    % Level 1: Average by exact Lat/Lon
                    [uni_locs_ok, ~, ic_ok] = unique([lat_ok, lon_ok], 'rows');
                    oman_uni_ok = accumarray(ic_ok, oman_ok, [], @(x) mean(x, 'omitnan'));
                    ombg_uni_ok = accumarray(ic_ok, ombg_ok, [], @(x) mean(x, 'omitnan'));

                    counts_afterQC(1,nn,mm,counter) = length(uni_locs_ok);

                    % Level 2: Average by Grid Index
                    ind_min_OK = knnsearch(Mdl, uni_locs_ok);
                    [grid_idx_ok, ~, ic_grid_ok] = unique(ind_min_OK);
                    oman_grid_ok = accumarray(ic_grid_ok, oman_uni_ok, [], @(x) mean(x, 'omitnan'));
                    ombg_grid_ok = accumarray(ic_grid_ok, ombg_uni_ok, [], @(x) mean(x, 'omitnan'));

                    oman_afterQC_vector_FULL(grid_idx_ok, nn, mm, counter) = oman_grid_ok;
                    ombg_afterQC_vector_FULL(grid_idx_ok, nn, mm, counter) = ombg_grid_ok;
                end
            end
        end
    end
    current_date = current_date + output_step;
end
fprintf('Data reading completed in %.2f seconds.\n', toc);

%% Plotting Time-Averaged Results
fprintf('Generating maps...\n');

abs_oman_beforeQC_vector_FULL = abs(oman_beforeQC_vector_FULL);
abs_ombg_beforeQC_vector_FULL = abs(ombg_beforeQC_vector_FULL);
abs_oman_afterQC_vector_FULL  = abs(oman_afterQC_vector_FULL);
abs_ombg_afterQC_vector_FULL  = abs(ombg_afterQC_vector_FULL);

oman_before_avg = mean(abs_oman_beforeQC_vector_FULL, 4, 'omitnan');
ombg_before_avg = mean(abs_ombg_beforeQC_vector_FULL, 4, 'omitnan');
oman_after_avg  = mean(abs_oman_afterQC_vector_FULL, 4, 'omitnan');
ombg_after_avg  = mean(abs_ombg_afterQC_vector_FULL, 4, 'omitnan');

for nn = 1:length(diag_source)
   switch diag_source{nn}
       case 'ims_snow'
          latlim = [0 90]; lonlim = [-180 180];
       case 'snocvr_snomad'
          latlim = [24 72]; lonlim = [-170 45];
       case 'sfcsno'
          latlim = [0 90]; lonlim = [-180 180];
       otherwise
          latlim = [-90 90]; lonlim = [-180 180];
   end

   for mm = 1:num_experiments
        fig = figure('Name', sprintf('%s_%s', exp_names{mm}, diag_source{nn}), 'Color', 'w');
        set(fig, 'Position', [100 100 1200 800]);
        t = tiledlayout(2,2, 'TileSpacing', 'compact', 'Padding', 'compact');
        
        data_to_plot = {ombg_before_avg(:,nn,mm), oman_before_avg(:,nn,mm), ...
                        ombg_after_avg(:,nn,mm), oman_after_avg(:,nn,mm)};
        titles = {'abs(OMBG) (Before JEDI QC)', 'abs(OMAN) (Before JEDI QC)', ...
                  'abs(OMBG) (After JEDI QC)', 'abs(OMAN) (After JEDI QC)'};
        
        for i = 1:4
            nexttile;
            worldmap(latlim, lonlim);
            current_data = data_to_plot{i};
            valid = ~isnan(current_data);
            domain_mean = mean(current_data, 'omitnan');

            if any(valid)
                scatterm(lat1D(valid), lon1D(valid), 4, current_data(valid), 'filled');
            end
            
            % Modular map finalization using cached geoshow boundaries
            fct_finalize_map(custom_colors, [cmin, cmax], '', afs, mask_dir);
            title(sprintf('%s: m=%.4f', titles{i}, domain_mean), 'Interpreter', 'none');
        end
        
        sgtitle(t, sprintf('Grid-averaged and Time-Averaged [%s]: %s | Exp: %s', date_range_str, diag_source{nn}, exp_names{mm}), ...
                'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'none');

        fname = fullfile(output_dir, sprintf('avg_diag_%s_%s.jpg', diag_source{nn}, exp_names{mm}));
        saveas(fig, fname);
        fprintf('Saved: %s\n', fname);
    end
end

%% Plotting Observation Counts Time Series
fprintf('Generating station count time series...\n');

for nn = 1:length(diag_source)
    for mm = 1:num_experiments
        fig = figure('Name', sprintf('Counts_%s_%s', exp_names{mm}, diag_source{nn}), 'Color', 'w');
        set(fig, 'Position', [100 100 1000 500]);
       
        c_before = squeeze(counts_beforeQC(1,nn,mm,:));
        c_after  = squeeze(counts_afterQC(1,nn,mm,:));

        plot(date_vector, c_before, '-o', 'LineWidth', 1.5, 'MarkerSize', 4);
        hold on;
        plot(date_vector, c_after, '-s', 'LineWidth', 1.5, 'MarkerSize', 4);
        hold off;
        ylabel('Total Stations');
        
        grid on;
        datetick('x', 'mm/dd HH', 'keepticks');
        xlabel('Date/Time (MM/DD HH)');
        title(sprintf('Station Counts (by unique lat/lon pairs): %s | Exp: %s', diag_source{nn}, exp_names{mm}), 'Interpreter', 'none');
        legend({'Before JEDI QC', 'After JEDI QC'}, 'Location', 'best');
        
        fname = fullfile(output_dir, sprintf('station_counts_timeseries_%s_%s.jpg', diag_source{nn}, exp_names{mm}));
        saveas(fig, fname);
        fprintf('Saved: %s\n', fname);
    end
end

fprintf('All tasks finished successfully!\n');

%% =========================================================================
%% HELPER FUNCTIONS
%% =========================================================================

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
    if ~isempty(colorbar_string)
        set(get(hl, 'ylabel'), 'String', colorbar_string, 'FontSize', afs);
    end

    setm(gca, 'LabelFormat', 'signed', 'FontSize', afs);
    if ~isempty(cached_land_lat), geoshow(cached_land_lat, cached_land_lon, 'Color', 'black'); end
    if ~isempty(cached_intl_lat), geoshow(cached_intl_lat, cached_intl_lon, 'Color', 'black', 'LineWidth', 0.5); end
    if ~isempty(cached_usa_lat),  geoshow(cached_usa_lat,  cached_usa_lon,  'Color', 'black', 'LineWidth', 0.5); end
end

function [lat, lon] = load_shp_coords(mask_dir, filename)
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
