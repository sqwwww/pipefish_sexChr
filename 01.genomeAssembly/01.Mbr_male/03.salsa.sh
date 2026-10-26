#!/bin/bash

#SBATCH --job-name=HJ_salsa
#SBATCH --partition=cpu
#SBATCH --ntasks-per-node=12
#SBATCH --mem=100g
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH --time=330:00:00


ref=HJ_hifi.fasta
fq1=fastq/HJ-M_R1.fastq.gz
fq2=fastq/HJ-M_R2.fastq.gz

module load bedtools-2.29.2
module load bwa-0.7.17
module load samtools-1.11


### prepare data
##samtools faidx $ref
##bwa index $ref
#bwa mem -t 12 -M $ref $fq1 $fq2 | samtools view -Sbh -q 20 > alignment.bam
#bamToBed -i alignment.bam > alignment.bed
#sort -k 4 alignment.bed > tmp && mv tmp alignment.bed
#
## run salsa
source /share/home/zhou3lab/suqianwei/miniconda3/bin/activate /share/home/zhou3lab/suqianwei/miniconda3/envs/mamba/envs/salsa2_reqs
/data/home/zhou3lab/suqianwei/tools2/SALSA-2.3/bin/run_pipeline.py --assembly $ref --length ${ref}.fai \
--bed alignment.bed --enzyme GATC --output out

# quast
source /share/home/zhou3lab/suqianwei/miniconda3/bin/activate quast
quast.py -t 10 out/scaffolds_FINAL.fasta
