#!/bin/bash

#SBATCH --job-name=era5land_regrid
#SBATCH -t 00:30:00
#SBATCH -A da-cpu
#SBATCH --qos=batch
#SBATCH -o era5_regrid.out
#SBATCH -e era5_regrid.err
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1 

source ../sorc/land_mods
./run_era5land_snow.sh



