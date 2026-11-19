#!/bin/bash
#SBATCH --job-name=pca
#SBATCH --time=70:00:00
#SBATCH --partition=cpu
#SBATCH -N 1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=10G
#SBATCH --output=%j.out
#SBATCH --error=%j.err

# SNP-based principal component analysis (Supplementary Figure S10).
# This script was run separately for four SNP sets:
# (A) genome-wide biallelic SNPs;
# (B) biallelic SNPs on Chr1;
# (C) male-specific SNPs on Chr1;
# (D) biallelic SNPs on Chr2.
#
# For each SNP set, prepare snp.sample and snp.haps in a separate
# working directory and run this script from that directory.
# Requirements: Perl, PLINK, and getsnpplink.pl.

# Convert SAMPLE/HAPS files to PLINK MAP/PED format.
# Place getsnpplink.pl in the working directory.
perl getsnpplink.pl snp.sample snp.haps snp.map snp.ped

# Convert MAP/PED files to binary PLINK format.
plink --file snp --out snp --make-bed --allow-extra-chr --allow-no-sex

# Export TPED/TFAM files (not required for the PCA step below).
plink --bfile snp --autosome-num 19 --allow-no-sex \
    --recode 12 transpose --output-missing-genotype 0 --out snp

# Perform PCA.
# Outputs: snp.eigenvec (PC coordinates) and snp.eigenval (eigenvalues).
plink --bfile snp --autosome-num 19 --allow-no-sex --pca --out snp