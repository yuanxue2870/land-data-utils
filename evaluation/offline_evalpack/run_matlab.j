#!/bin/bash -le
#SBATCH --job-name=matlab
#SBATCH --account=da-cpu
#SBATCH --qos=batch
#SBATCH --nodes=1
#SBATCH --mem=300000
#SBATCH -t 04:30:00

module load matlab
# Uncomment lines below to run any experiment(s) needed
#matlab -nodisplay -nosplash -nodesktop -r "run('Evaluate_snow_depth.m');exit;"
#matlab -nodisplay -nosplash -nodesktop -r "run('Evaluate_scf.m');exit;"
#matlab -nodisplay -nosplash -nodesktop -r "run('Analyze_diag.m');exit;"
#matlab -nodisplay -nosplash -nodesktop -r "run('Analyze_incr.m');exit;"
