#!/bin/bash

#BATCH --job-name=fv3_mapping_land
#SBATCH -t 07:55:00
#SBATCH -A bigmem
#SBATCH -A da-cpu
#SBATCH --qos=batch
#SBATCH -o mapping_snow.out
#SBATCH -e mapping_snow.out
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1 

source mods_bash

./create_land_mapping.exe



