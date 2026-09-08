#!/bin/bash

#SBATCH --job-name=land_mapping
#SBATCH -t 07:30:00
#SBATCH -A da-cpu
#SBATCH --qos=batch
#SBATCH -o mapping_snodas.out
#SBATCH -e mapping_snodas.err
#SBATCH --nodes=2
#SBATCH --tasks-per-node=40

export OMP_NUM_THREADS=80
source ../land_mods
./create_land_mapping.exe



