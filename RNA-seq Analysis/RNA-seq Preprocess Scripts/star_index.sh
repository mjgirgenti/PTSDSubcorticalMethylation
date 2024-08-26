#!/bin/bash

#SBATCH --job-name=create-genome-index-human
#SBATCH --mem=60G
#SBATCH -c 8
#SBATCH --mail-type=ALL

PATH_TO_FASTA=/home/tpn5/gibbs/HPC/RefGenome/Homo_sapiens.GRCh38.dna.primary_assembly.fa # Path to genome assembly
PATH_TO_GTF=/home/tpn5/gibbs/HPC/RefGenome/Homo_sapiens.GRCh38.96.gtf

OUT_FOLDER=/home/tpn5/gibbs/HPC/RefGenome/STARindices

module load STAR

STAR --runThreadN 8 --runMode genomeGenerate --genomeDir $OUT_FOLDER --genomeFastaFiles $PATH_TO_FASTA --sjdbGTFfile $PATH_TO_GTF
