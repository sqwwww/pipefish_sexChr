#!/bin/bash
#SBATCH --job-name=fst
#SBATCH --time=12:00:00 # hh:mm:ss
#SBATCH --partition=short
#SBATCH -N 1
#SBATCH --ntasks-per-node=2
#SBATCH --mem=70G
#SBATCH --output=%j.out
#SBATCH --error=%j.err



vcftools --gzvcf ref.snp.final.vcf.gz --weir-fst-pop population1.txt --weir-fst-pop population2.txt --out Fst_bin --fst-window-size 100000 --fst-window-step 50000 
