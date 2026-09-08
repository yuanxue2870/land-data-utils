%% config_inputs.m
%% Note: Gemini and Github Copilot were used to assist with developing this code. The code has been reviewed, edited, and validated by NWS staff.

% Centralized configuration for snow depth and SCF evaluation and analysis scripts

% =========================================================================
% 1. CORE COMMON INPUTS (Shared by all scripts)
% =========================================================================
res = 1152; %tile/land/model grid resolution
ocn_res = 'mx025'; %ocean grid resolution

% User's offline workflow output directory
directory = '/scratch3/NCEPDEV/stmp/Yuan.Xue/eval_test/';

% Eval start date and end date
start_date = datenum(2024, 11, 1);
finish_date = datenum(2025, 2, 1);

% =========================================================================
% 2. SCRIPT-SPECIFIC CONFIGURATIONS
% =========================================================================

% --- Settings for Evaluate_snow_depth(scf).m ---
eval_output_step = 1; %unit: day
eval_exp_names = {'ol'; 'ori'; 'rejectA'; 'rejectB'; 'rejectC'}; %Note: put the baseline as the first experiment in exp_names
threshold_num_days = 30; %minimum number of samples needed to generate goodness-of-fit metrics

% --- Settings for Analyze_diag.m and Analyze_incr.m ---
diag_output_step = 0.25; %unit: day; 0.25 means every 6 hours
incr_exp_names = {'ori'; 'rejectA'; 'rejectB'; 'rejectC'};
diag_source = {'snocvr_snomad'}; %options to choose(based on user's DA experiments): ims_snow, snocvr_snomad, sfcsno
region_of_interest = 'US_plus_Nordic';%options to choose: ND, US, CONUS, US_plus_Nordic
