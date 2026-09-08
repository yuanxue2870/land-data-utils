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
- **Land fraction mask files** (pre-computed .mat files with grid information, currently can be accessed from Yuan Xue's folder)

### Step 1: Configure Input Parameters

Edit `config_inputs.m` to set up your evaluation:

### Step 2: Run Analysis Scripts

Each script can be run independently and/or via the job scheduler:

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

## Support & Questions

For questions or issues:
1. Review diagnostic output in MATLAB console
2. Contact the package maintainer (Yuan Xue)

---

**Last Updated**: May 11, 2026
**Version**: 1.0
**Status**: Production Ready
