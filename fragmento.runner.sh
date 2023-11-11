#!/bin/bash
# the following are sbatch parameters
#SBATCH -t 72:00:00
#SBATCH -p himem
#SBATCH --mem=60G
#SBATCH -c 1
#SBATCH -N 1
#SBATCH -o %x-%j.out

module load R/4.1

echo 'Running fragment analysis script'
Rscript fragmento.custom.R
echo 'Finished running fragment analysis script!'
