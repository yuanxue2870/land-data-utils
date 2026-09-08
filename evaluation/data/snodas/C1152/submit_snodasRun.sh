#!/bin/bash

#BATCH --job-name=snodas_run
#SBATCH -t 03:30:00
#SBATCH -A da-cpu
#SBATCH --qos=batch
#SBATCH -o snodas_run.out
#SBATCH -e snodas_run.err
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1 

source ../sorc/land_mods
./run_snodas.sh
