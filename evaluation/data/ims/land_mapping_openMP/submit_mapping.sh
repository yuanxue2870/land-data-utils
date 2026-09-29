#!/bin/bash

#BATCH --job-name=fv3_mapping_land
#SBATCH -t 07:55:00
#SBATCH -A da-cpu
#SBATCH --qos=batch
#SBATCH -o mapping_snow.out
#SBATCH -e mapping_snow.err
#SBATCH --nodes=2
#SBATCH --tasks-per-node=40

export OMP_NUM_THREADS=80

./create_land_mapping.exe



