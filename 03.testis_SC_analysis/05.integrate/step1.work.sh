#!/bin/bash
#SBATCH -J sc
#SBATCH -p short
#SBATCH -N 1
#SBATCH --mem=40G
#SBATCH --ntasks-per-node=4
#SBATCH --output=%J.out
#SBATCH --error=%J.err


source /share/home/zhou3lab/suqianwei/miniconda3/bin/activate /share/home/zhou3lab/suqianwei/miniconda3/envs/mamba/envs/seurat5
#Rscript combinedCintHa.R
#Rscript combinedHa_zbf_Cint.R
#Rscript combinedMbrDre.R 

#Rscript combinedHa_zbf_Cint_Mbr.R
Rscript combinedSeurat.R
