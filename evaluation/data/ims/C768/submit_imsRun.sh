#!/bin/bash

#BATCH --job-name=ims_run
#SBATCH -t 07:55:00
#SBATCH -A bigmem
#SBATCH -A da-cpu
#SBATCH --qos=batch
#SBATCH -o ims_run.out
#SBATCH -e ims_run.out
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1 

./run_ims_snow.sh



