#!/bin/bash

#SBATCH --job-name=fv3_mapping
#SBATCH -t 01:30:00
#SBATCH -A da-cpu
#SBATCH --qos=batch
#SBATCH -o mapping_snodas.out
#SBATCH -e mapping_snodas.err
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1

source ../land_mods
./create_fv3_mapping.exe



