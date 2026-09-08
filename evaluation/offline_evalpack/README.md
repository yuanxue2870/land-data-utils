# Offline Evaluation Package (offline_evalpack)

A comprehensive MATLAB-based evaluation and analysis toolkit for offline land data assimilation (DA) experiments. This package provides tools to evaluate snow depth and snow cover fraction (SCF) model performance, analyze JEDI diagnostic files, and compare data assimilation increments across multiple experiments.

## Overview

The offline evaluation package is designed for the land data assimilation workflow and includes:

- **Snow Depth & SCF Evaluation**: Compare model performance against observations
- **Diagnostic Analysis**: Analyze observation minus forecast (O-F) statistics from JEDI
- **Increment Analysis**: Visualize and compare DA increments across experiments
- **Automated Processing**: Batch processing capabilities via job scheduler integration

## Directory Structure

```
offline_evalpack/
├── README.md                    # This file
├── config_inputs.m              # Central configuration file (user inputs)
├── Evaluate_snow_depth.m        # Snow depth evaluation script
├── Evaluate_scf.m               # Snow cover fraction (SCF) evaluation script
├── Analyze_diag.m               # Diagnostic file analysis
├── Analyze_incr.m               # Data assimilation increment analysis
└── run_matlab.j                 # SLURM job submission script
```

## Getting Started

### Prerequisites

- **MATLAB R2019a or later** (for compatibility with mapping toolbox and modern syntax)
- **MATLAB Mapping Toolbox** (required for map visualizations and shapefile I/O)
- **NetCDF files** from offline land DA workflow (diagnostic and increment files)
- **Land fraction mask files** (pre-computed .mat files with grid information)

### Step 1: Configure Input Parameters

Edit `config_inputs.m` to set up your evaluation:

```matlab
%% Grid Configuration
res = 1152;              % Tile/model grid resolution (384, 768, or 1152)
ocn_res = 'mx025';       % Ocean grid resolution (currently only mx025 supported)

%% Data Paths
directory = '/path/to/your/offline/workflow/output/';

%% Time Period
start_date = datenum(2024, 11, 1);    % Start date (YYYY, MM, DD)
finish_date = datenum(2025, 2, 1);    % End date (YYYY, MM, DD)

%% Evaluation Settings
eval_output_step = 1;                          % Daily evaluation output
eval_exp_names = {'ol'; 'ori'; 'rejectA'};   % Baseline experiment first
threshold_num_days = 30;                       % Min samples for statistics

%% Diagnostic Analysis Settings
diag_output_step = 0.25;                       % 6-hourly analysis
diag_source = {'snocvr_snomad'};              % Observation data source
region_of_interest = 'US_plus_Nordic';        % Spatial domain

%% Increment Analysis Settings
incr_exp_names = {'ori'; 'rejectA'; 'rejectB'};  % Experiments to compare
```

**Key Configuration Notes:**
- **Baseline experiment**: First entry in `eval_exp_names` is treated as the reference/baseline
- **Grid resolution**: Must be one of `[384, 768, 1152]`
- **Time format**: Use MATLAB `datenum()` (days since Jan 1, 0000)
- **Output step**: Specify in days (0.25 = 6 hours, 1 = daily)
- **Observation sources**: Options include `'ims_snow'`, `'snocvr_snomad'`, `'sfcsno'`
- **Regions**: `'ND'`, `'US'`, `'CONUS'`, `'US_plus_Nordic'`

### Step 2: Run Analysis Scripts

Each script can be run independently or via the job scheduler:

#### Interactive Mode (single experiment)

```matlab
% In MATLAB command window:
run('Analyze_diag.m')
run('Analyze_incr.m')
run('Evaluate_snow_depth.m')
```

#### Batch Mode (via SLURM scheduler)

Edit `run_matlab.j` to uncomment desired scripts:

```bash
# Uncomment the analysis you want to run:
matlab -nodisplay -nosplash -nodesktop -r "run('Evaluate_snow_depth.m');exit;"
```

Then submit:
```bash
sbatch run_matlab.j
```

## Script Descriptions

### `config_inputs.m`
Central configuration file containing all user-adjustable parameters. This file is sourced by all analysis scripts to ensure consistency across the workflow.

**Key outputs**: Defines variables used by all other scripts
- `res`, `ocn_res`: Grid resolutions
- `directory`: Workflow output directory
- `start_date`, `finish_date`: Evaluation time range
- `exp_names`, `eval_exp_names`: Experiment identifiers

### `Analyze_diag.m`
Analyzes JEDI diagnostic files containing observation-minus-forecast (O-F) statistics.

**Inputs**:
- NetCDF diagnostic files: `diag_*.nc`
- Mask file: `C{res}_hr3_{ocn_res}_era5_input.mat`

**Outputs**:
- `Images_diag_analysis/avg_diag_*.jpg`: Time-averaged O-F maps
- `Images_diag_analysis/station_counts_timeseries_*.jpg`: Observation count time series

**Key functions**:
- Double averaging: Unique lat/lon averaging + grid-cell averaging
- QC filtering: Applies JEDI QC flags before/after quality control
- Visualization: Side-by-side maps of OMBG and OMAN (before and after JEDI QC)

### `Analyze_incr.m`
Analyzes data assimilation increments and creates animated comparisons.

**Inputs**:
- NetCDF increment files: `*.snow_increment.sfc_data.tile*.nc`
- Mask file with 2D lat/lon coordinates

**Outputs**:
- `snodl_increment_comparison_{roi}.gif`: Animated increment comparison across experiments
- Console plots with min/max increment markers

**Key functions**:
- Tile-based processing (e.g., 14 tiles for C1152 resolution)
- Regional cropping and visualization
- Blue-to-Red colormap (symmetric around zero)

### `Evaluate_snow_depth.m`
Comprehensive snow depth evaluation against observations.

**Inputs**:
- Model output and observation data files
- Mask file

**Outputs**:
- Goodness-of-fit metrics (bias, RMSE, correlation, etc.)
- Evaluation maps and statistics
- Baseline-relative comparisons

### `Evaluate_scf.m`
Comprehensive snow cover fraction (SCF) evaluation against observations.

**Inputs**:
- Model output and observation data files
- Mask file

**Outputs**:
- SCF evaluation metrics
- Contingency tables and error rates
- Performance comparison plots

### `run_matlab.j`
SLURM batch job submission script for high-performance computing clusters.

**Configuration**:
- Account: `da-cpu`
- QoS: `batch`
- Memory: 300 GB
- Time limit: 4.5 hours
- Nodes: 1

To use: Edit the `matlab` command lines to uncomment desired scripts, then run:
```bash
sbatch run_matlab.j
```

## Expected Output Structure

```
offline_evalpack/
├── Images_diag_analysis/
│   ├── avg_diag_snocvr_snomad_ori.jpg
│   ├── avg_diag_snocvr_snomad_rejectA.jpg
│   ├── station_counts_timeseries_snocvr_snomad_ori.jpg
│   ├── station_counts_timeseries_snocvr_snomad_rejectA.jpg
│   └── ...
├── snodl_increment_comparison_US_plus_Nordic.gif
├── evaluation_results/
│   └── [evaluation metrics and plots]
└── config_inputs.m
```

## Data Requirements

### Offline Workflow Directory Structure

The evaluation scripts expect output from your offline DA workflow organized as:

```
{directory}/
├── ol/
│   ├── DA/jedi_anl/
│   │   ├── diag_*.nc
│   │   └── *.snow_increment.sfc_data.tile*.nc
│   └── [other output files]
├── ori/
│   ├── DA/jedi_anl/
│   │   └── [same structure as above]
│   └── ...
├── rejectA/
└── rejectB/
```

### Mask Files

Mask files must be located in the directory specified by `mask_dir`:

```
{mask_dir}/
├── C384_hr3_mx025_era5_input.mat
├── C768_hr3_mx025_era5_input.mat
├── C1152_hr3_mx025_era5_input.mat
└── [shapefile data for map overlays]
```

**Mask file contents**:
- `C{res}_lat1D`, `C{res}_lon1D`: 1D coordinate arrays
- `C{res}_lat2D_FULL`, `C{res}_lon2D_FULL`: 2D coordinate arrays per tile
- Land fraction masks

## Common Workflow

### Quick Start Example

```matlab
% 1. Set configuration
edit config_inputs.m    % Adjust paths, dates, experiments

% 2. Run diagnostic analysis
run('Analyze_diag.m')   % Outputs: Images_diag_analysis/*.jpg

% 3. Run increment analysis
run('Analyze_incr.m')   % Outputs: snodl_increment_comparison_*.gif

% 4. Run full evaluations (optional, more computationally intensive)
run('Evaluate_snow_depth.m')
run('Evaluate_scf.m')
```

### Troubleshooting

**Error: "Mask file does not exist"**
- Check that `mask_dir` path in `config_inputs.m` is correct
- Verify mask files exist at the specified location

**Error: "Invalid resolution"**
- Ensure `res` is one of: 384, 768, or 1152
- Check that mask file matches the specified resolution

**Error: "No data found for date"**
- Verify `start_date` and `finish_date` are valid
- Confirm experiment directories exist in the offline workflow output path
- Check that diagnostic/increment files follow the expected naming convention

**Out of memory errors**
- Reduce time range (`finish_date - start_date`)
- Decrease `eval_output_step` (fewer dates)
- Reduce grid resolution if possible

**Missing shapefiles for maps**
- Download NOAA Natural Earth or USGS shapefiles
- Place in `mask_dir`
- Required files: `landareas.shp`, `intl_boundaries.shp`, `usastatelo.shp`

## Algorithm Notes

### Double Averaging (Analyze_diag.m)

1. **Level 1**: Average by unique (lat, lon) pairs (handles repeated observations)
2. **Level 2**: Map observations to nearest model grid cell using KD-tree
3. **Level 3**: Average values within grid cells

This approach ensures proper weighting of observations with multiple measurements at the same location.

### QC Filtering

- **Before QC**: All observations (EffectiveQC0 and EffectiveQC1 may not be applied)
- **After QC**: Only observations passing both QC0==0 and QC1==0 checks

### Blue-to-Red Colormap (Analyze_incr.m)

- **Blue**: Negative values (analysis reducing snow depth)
- **White**: Zero (no change)
- **Red**: Positive values (analysis increasing snow depth)
- Symmetric scaling around zero for balanced visualization

## Development Notes

**Code Quality**:
- Code assisted by Gemini and GitHub Copilot
- All code reviewed, edited, and validated by NWS staff
- Extensive error checking and validation

**Future Enhancements**:
- Uncertainty quantification
- Multi-variable evaluation (soil moisture, temperature, etc.)
- Real-time monitoring capabilities
- Database storage of evaluation metrics

## Citation & Attribution

This package was developed for the NOAA National Weather Service land data assimilation workflow. 

**Authors**:
- Yuan Xue (NOAA/NWS)
- Assisted by Gemini and GitHub Copilot

## License

See repository LICENSE file.

## Support & Questions

For questions or issues:
1. Check the troubleshooting section above
2. Review diagnostic output in MATLAB console
3. Contact the package maintainer (Yuan Xue)

---

**Last Updated**: May 11, 2026
**Version**: 1.0
**Status**: Production Ready
