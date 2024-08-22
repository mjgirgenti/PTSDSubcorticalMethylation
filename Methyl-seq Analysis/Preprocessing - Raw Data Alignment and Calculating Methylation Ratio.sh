
## fastp algorithm:
## Input files:
##	r1 = fastqdir + "/{sample}_R1.fastq.gz"
##      	r2 = fastqdir + "/{sample}_R2.fastq.gz" 
## Output files:
##	r1 = temp(fastpdir + "/{sample}_R1.fastq.gz"),
##      	r2 = temp(fastpdir + "/{sample}_R2.fastq.gz"),
##      	json_report = fastpdir + "/json_reports/{sample}.json",
##      	html_report = fastpdir + "/html_reports/{sample}.html"
## Versions: 
##	fastp ==0.20.0 
## Rules:
ml miniconda; conda activate fastp; 
fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} --json {output.json_report} --html {output.html_report} --thread 10 --detect_adapter_for_pe



## BSMAP algorithm:
## Input files:
##	r1 = fastpdir + "/{sample}_R1.fastq.gz"
##      	r2 = fastpdir + "/{sample}_R2.fastq.gz"
## Output files:
##      	bam = temp(bsmapdir + "/{sample}.bam")
##      	bai = temp(bsmapdir + "/{sample}.bai")
## Versions:
## 	bsmap ==2.90
##	python ==2.7.15
##	samtools ==0.1.19
## 	openssl=1.1.1w
## Rules:
ml miniconda; conda activate BSMAP; 
(set -x; bsmap -a {input.r1} -b {input.r2} -d Homo_sapiens.GRCh38.dna.primary_assembly.fa -p 20 -r 0 -w 2 -v 5 | samtools view -f2 -Shu - | samtools sort -l9 -@5 - {params.prefix} && samtools index {output.bam} {params.prefix}.bai) &> {log}

## The reference genome is hg38 version of human genome, Ensembl version 96 downloaded from: http://ftp.ensembl.org/pub/release-96/fasta/homo_sapiens/dna/



## MarkDuplicate algorithm:
## Input files:
## 	bsmapdir + "/{sample}.bam"
## Output files:
##      bam = temp(markdupdir + "/{sample}.bam")
##      bai = temp(markdupdir + "/{sample}.bai")
##      metrics = markdupdir + "/metrics/{sample}.metrics"
## Versions:
##	picard ==2.21.6
## Rules:
ml miniconda; conda activate picard; 
(set -x; picard MarkDuplicates -Xmx8G INPUT={input} OUTPUT={output.bam}METRICS_FILE={output.metrics} REMOVE_DUPLICATES=true CREATE_INDEX=true) &> {log}



## methratio algorithm in BSMAP (same environment with BSMAP algorithm):
## Input files:
##	bam = markdupdir + "/{sample}.bam"
## Output files:
##      	temp(methratiodir + "/{sample}-methratio.tsv")
## Rules: 
ml miniconda; conda activate BSMAP; 
(set -x; methratio.py -d Homo_sapiens.GRCh38.dna.primary_assembly.fa {input.bam} -o {output}) &> {log}

## The reference genome is hg38 version of human genome, Ensembl version 96 downloaded from: http://ftp.ensembl.org/pub/release-96/fasta/homo_sapiens/dna/


