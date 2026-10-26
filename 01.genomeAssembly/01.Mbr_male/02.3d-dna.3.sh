#!/bin/bash
#SBATCH --job-name=3d-dna
#SBATCH --partition=cpu
#SBATCH -N 1
#SBATCH --mem=300G
#SBATCH --ntasks-per-node=24
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH --time=330:00:00

source /share/home/zhou3lab/suqianwei/miniconda3/bin/activate 3d-dna
module load samtools-1.11
module load bwa-0.7.17

#bwa index -p reference/HJ_hifi_hap2.fasta reference/HJ_hifi_hap2.fasta


#python /data/home/zhou3lab/suqianwei/tools2/juicer/misc/generate_site_positions.py MboI HJ reference/HJ_hifi_hap2.fasta
#awk 'BEGIN{OFS="\t"}{print $1, $NF}' HJ_MboI.txt > HJ.chrom.sizes

#/data/home/zhou3lab/suqianwei/tools2/juicer/scripts/juicer.sh \
#-g HJ -q cpu -Q 330:00:00 -l cpu -L 330:00:00 \
#-D /data/home/zhou3lab/suqianwei/tools2/juicer \
#-z reference/HJ_hifi_hap2.fasta \
#-y HJ_MboI.txt \
#-p HJ.chrom.sizes \
#--assembly \
#-s MboI -t 24 -T 24

/data/home/zhou3lab/suqianwei/tools2/3d-dna/run-asm-pipeline.sh -r 0 reference/HJ_hifi_hap2.fasta aligned/merged_nodups.txt

source /share/home/zhou3lab/suqianwei/miniconda3/bin/activate quast
quast.py -t 10 HJ_hifi_hap2.FINAL.fasta
