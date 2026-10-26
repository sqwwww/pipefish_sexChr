#!/bin/bash
#SBATCH --job-name=kmer
#SBATCH --partition=short
#SBATCH --cpus-per-task=18
#SBATCH --mem=20G
#SBATCH --output=kmer.o
#SBATCH --error=kmer.e
#conda activate merqury
###illumina 15 17 21; hifi less than 25
yak count -o k31.yak -k 31 -b 37 -t 18 <(cat BJ_female_DNA_*.fq) <(zcat BJ_male_DNA_*.fq.gz)
yak count -o k21.yak -k 21 -b 37 -t 18 <(cat BJ_female_DNA_*.fq) <(zcat BJ_male_DNA_*.fq.gz)
