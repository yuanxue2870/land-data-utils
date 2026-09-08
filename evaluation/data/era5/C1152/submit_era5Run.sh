#!/bin/bash

#SBATCH --job-name=era5_regrid
#SBATCH -t 04:00:00
#SBATCH -A da-cpu
#SBATCH --qos=batch
#SBATCH -o era5_run.out
#SBATCH -e era5_run.err
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1 

./run_era5_snow.sh



