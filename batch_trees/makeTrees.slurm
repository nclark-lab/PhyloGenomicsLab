#!/bin/bash
#SBATCH --job-name=makeTrees            
#SBATCH -t 72:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --partition=htc            
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=nclark@pitt.edu
#SBATCH -o slurm-%j-%N.o
#SBATCH -e slurm-%j-%N.e

module purge
module load r/4.5.0
cd /ihome/nclark/nclark/forNathan
Rscript makeTrees.R
