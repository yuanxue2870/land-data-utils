#!/bin/bash

#BATCH --job-name=ims_run
#SBATCH -t 07:55:00
#SBATCH -A da-cpu
#SBATCH --qos=batch
#SBATCH -o ims_run.out
#SBATCH -e ims_run.err
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1 

source ../sorc/land_mods
./run_ims_snow.sh



