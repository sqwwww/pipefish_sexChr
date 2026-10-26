#!/bin/bash
#SBATCH --job-name=convert
#SBATCH --time=12:00:00 # hh:mm:ss
#SBATCH --partition=short
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=10G
#SBATCH --output=%j.out
#SBATCH --error=%j.err


plink2 \
  --export vcf \
  --haps snp.haps \
  --out snp.haps.out

