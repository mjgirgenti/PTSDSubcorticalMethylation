#!/bin/bash

#SBATCH -p day -t 1-0
#SBATCH --job-name=check-quality-AMG-HPC
#SBATCH -c 5
#SBATCH --mail-type=ALL

PATH_TO_FASTQ=/home/tpn5/lab/methylation_shared_folder/RNAseq/RawFastq
OUTDIR=/home/tpn5/lab/methylation_shared_folder/RNAseq/Results

#-------------------

module load FastQC
module load MultiQC

fastqc -t 5 $PATH_TO_FASTQ/*/Unaligned/*.fastq.gz

multiqc $PATH_TO_FASTQ -o $OUTDIR
