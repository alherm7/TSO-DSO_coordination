#!/bin/sh
#BSUB -J PCC_optimzer
#BSUB -q elektro
#BSUB -n 1
#BSUB -R "rusage[mem=2GB]"
#BSUB -M 10GB
#BSUB -W 20:00
#BSUB -u alherm@dtu.dk
#BSUB -B
#BSUB -N
#BSUB -o output_alex_run1.out
#BSUB -e error_alex_run1.err
#BSUB -R "span[hosts=1]"

module load cvx
module load mosek/9.2
##module load gurobi/8.1.1

matlab -nodisplay -r RUN_PCC_optim -logfile PCC_optim_logfile_output


##BSUB -m "n-62-21-94"
