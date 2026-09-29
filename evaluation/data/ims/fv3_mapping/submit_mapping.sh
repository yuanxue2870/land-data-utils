#!/bin/bash

#BATCH --job-name=fv3_mapping
#SBATCH -t 08:00:00
#SBATCH -A da-cpu
#SBATCH --qos=batch
#SBATCH -o mapping.out
#SBATCH -e mapping.err
#SBATCH --nodes=15
#SBATCH --tasks-per-node=40

source mods_bash 

./create_fv3_mapping.exe



