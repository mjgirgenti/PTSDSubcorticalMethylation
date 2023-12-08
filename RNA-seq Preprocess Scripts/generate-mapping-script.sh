#!/bin/bash

PATH_TO_FASTQ=/home/tpn5/lab/methylation_shared_folder/RNAseq/RawFastq
PATH_TO_GEN_IDX=/home/tpn5/gibbs/HPC/RefGenome/STARindices
OUTDIR=/home/tpn5/lab/methylation_shared_folder/RNAseq/MapFromRaw

PARALLEL='True' # True if we want to use dSQ for parallel mapping

#------------------------------------------------

if [[ -f "star_map.sh" ]]; then # Remove if the file already exist
rm star_map.sh 
fi

touch star_map.sh

if [ $PARALLEL == False ]; then
echo "#!/bin/bash" >> star_map.sh
echo "#SBATCH -p day -t 1-0 --job-name mapping-AMG/HPC --mem 30G -c 10 --mail-type=ALL" >> star_map.sh
fi

for sample in $PATH_TO_FASTQ/*/Unaligned; do 
inpath1=$(readlink -f $sample/*R1_001.fastq.gz | paste -sd, -) 
inpath2=$(readlink -f $sample/*R2_001.fastq.gz | paste -sd, -) 
outname=$(basename $(dirname $sample))
printf '%s\n' "module load STAR ; STAR --runThreadN 10 --runMode alignReads --genomeDir $PATH_TO_GEN_IDX --readFilesIn $inpath1 $inpath2 --readFilesCommand gunzip -c --outFileNamePrefix $OUTDIR/$outname --outSAMtype BAM SortedByCoordinate" >> star_map.sh 
done

if [ $PARALLEL == True ]; then
module load dSQ
dsq --job-file star_map.sh --job-name mapping-AMG/HPC -c 10 --mail-type=ALL
fi
