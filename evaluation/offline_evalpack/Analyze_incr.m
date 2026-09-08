%% The script is used to analyze snodl increments
%% Author: Yuan Xue, Gemini, Github Copilot
%% Date: 05/11/2026
%% Note: Gemini and Github Copilot were used to assist with developing this code. The code has been reviewed, edited, and validated by NWS staff.

clc; clear;

%% Constants and fixed eval directories
LARGE_SNOW_LIMIT = 1e12;
afs = 13;
cmin = -15; % lowest incr in plots
cmax = 15;  % highest incr in plots
custom_colors_diff = b2r(cmin, cmax);

%% User inputs
config_inputs; % Runs the central config and loads all variables

% Map script-specific variations
output_step = diag_output_step;
exp_names   = incr_exp_names;

%% Derived variables
if ischar(exp_names)
    exp_names = strtrim(strsplit(exp_names, ','));
elseif isstring(exp_names)
    exp_names = cellstr(strtrim(exp_names));
end

num_experiments = numel(exp_names);
total_steps = numel(start_date:output_step:finish_date);
date_vector = start_date:output_step:finish_date;

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
if isnan(output_step) || output_step <= 0
    error('Output step must be a positive number.');
end
if ~exist(directory, 'dir')
    error('Experiment directory does not exist: %s', directory);
end
if isempty(exp_names) || ~iscellstr(exp_names)
    error('Experiment names must be a non-empty cell array of strings.');
end

%% Load land fraction mask
mask_dir  = '/scratch3/NCEPDEV/stmp/Yuan.Xue/stats/';
mask_file = sprintf('%sC%d_hr3_%s_era5_input.mat', mask_dir, res, ocn_res);
if ~exist(mask_file, 'file')
    error('Mask file does not exist: %s', mask_file);
end
maskData = load(mask_file);
lat2D = maskData.(sprintf('C%d_lat2D_FULL', res));
lon2D = maskData.(sprintf('C%d_lon2D_FULL', res));
num_tiles = size(lat2D, 3);

% Define ROI Bounding Box
switch region_of_interest
   case 'ND'
      lat_min = 45; lat_max = 49.5;
      lon_min = -105; lon_max = -96;
   case 'US'
      lat_min = 24; lat_max = 72;
      lon_min = -170; lon_max = -66;
   case 'CONUS'
      lat_min = 24; lat_max = 49.5;
      lon_min = -125; lon_max = -66;
   case 'US_plus_Nordic'
      lat_min = 24; lat_max = 72;
      lon_min = -170; lon_max = 45;
   otherwise
      error('The lat/lon box has not defined for region %s', region_of_interest);
end

fprintf('Starting run: res=%d, roi=%s, total_steps=%d, num_tiles=%d, experiments=%d\n', ...
        res, region_of_interest, total_steps, num_tiles, num_experiments);

%% GIF Configuration
gif_filename = ['snodl_increment_comparison_', region_of_interest, '.gif'];
delay_time = 0.6; % Seconds between frames

%% Process each date
tic;  % Start timing for data reading
current_date = start_date;
counter = 0;

% Set up a wide figure panel for side-by-side comparison
figure('Position', [100, 100, 1200, 500]);

while current_date <= finish_date
    counter = counter + 1;

    fprintf('Processing date: %s\n', datestr(current_date));
    yyyys = datestr(current_date, 'yyyy');
    mms   = datestr(current_date, 'mm');
    dds   = datestr(current_date, 'dd');
    hhs   = datestr(current_date, 'HH');

    clf; % Clear figure to prepare for new subplots
   
    %% Assemble experiment increments and plot side-by-side
    for mm = 1:num_experiments
        current_ex = exp_names{mm};
        incr_FULL = NaN(res, res, num_tiles);

        for i = 1:num_tiles
            file = sprintf('%s%s/DA/jedi_anl/%s%s%s.%s0000.snow_increment.sfc_data.tile%d.nc', ...
                          directory, current_ex, yyyys, mms, dds, hhs, i);
            if exist(file, 'file')
                incr = ncread(file, 'snodl');
                incr(incr > LARGE_SNOW_LIMIT) = NaN;
                incr_FULL(:,:,i) = incr;
            end
        end

        % Find min/max values and linear indices (ignoring NaNs)
        [min_val, min_idx] = min(incr_FULL, [], 'all', 'omitnan');
        [max_val, max_idx] = max(incr_FULL, [], 'all', 'omitnan');

        % Convert linear indices to (row, col, tile) subscripts
        [min_row, min_col, min_tile] = ind2sub(size(incr_FULL), min_idx);
        [max_row, max_col, max_tile] = ind2sub(size(incr_FULL), max_idx);

        subplot(1, num_experiments, mm);
        if strcmp(region_of_interest, 'US_plus_Nordic')
            ax = axesm('miller', 'MapLatLimit', [lat_min lat_max], 'MapLonLimit', [lon_min lon_max]);
        else
            ax = axesm('lambert', 'MapLatLimit', [lat_min lat_max], 'MapLonLimit', [lon_min lon_max]);
        end
        gridm on; framem on; tightmap; hold on;
        
        for i = 1:num_tiles
             pcolorm(squeeze(lat2D(:,:,i)), squeeze(lon2D(:,:,i)), squeeze(incr_FULL(:,:,i)));
        end

        h_min = plotm(lat2D(min_row,min_col,min_tile), lon2D(min_row,min_col,min_tile), ...
                      'kx', 'MarkerSize', 7, 'LineWidth', 1.5);
        h_max = plotm(lat2D(max_row,max_col,max_tile), lon2D(max_row,max_col,max_tile), ...
                      'k^', 'MarkerSize', 7, 'LineWidth', 1.5);

        % Centralized colormap, colorbar, clim, and cached geoshow overlays
        fct_finalize_map(custom_colors_diff, [cmin, cmax], '', afs, mask_dir);

        % Title and legend
        title(sprintf('Exp: %s\\newline%s\\newlineMin: %.4f, Max: %.4f', ...
                      current_ex, datestr(current_date,'yyyy-mm-dd HH:MM:SS'), ...
                      min_val, max_val), 'FontSize', afs);
        legend([h_min, h_max], {'Min Inc', 'Max Inc'}, ...
               'Location', 'southoutside', 'FontSize', afs-3, 'Orientation', 'horizontal');
    end

    drawnow;

    % --- GIF Conversion ---
    frame = getframe(gcf);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);

    if counter == 1
        imwrite(imind, cm, gif_filename, 'gif', 'Loopcount', inf, 'DelayTime', delay_time);
    else
        imwrite(imind, cm, gif_filename, 'gif', 'WriteMode', 'append', 'DelayTime', delay_time);
    end

    current_date = current_date + output_step;
end

processing_time = toc;
fprintf('Data reading completed in %.2f seconds.\n', processing_time);
fprintf('All tasks finished successfully!\n');

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
