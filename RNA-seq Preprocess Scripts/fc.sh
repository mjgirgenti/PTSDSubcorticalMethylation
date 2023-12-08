#!/bin/bash

#SBATCH -p day -t 1-0 -c 50
#SBATCH --job-name=feature-count
#SBATCH --mail-type=ALL

module load R/4.2.0-foss-2020b

Rscript fc.R
