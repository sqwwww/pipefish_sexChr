#!/bin/bash
#SBATCH -J velocyto
#SBATCH -p cpu
#SBATCH -N 1
#SBATCH --mem=30G
#SBATCH --ntasks-per-node=6
#SBATCH --output=%J.out
#SBATCH --error=%J.err




source /share/home/zhou3lab/suqianwei/miniconda3/bin/activate ~/miniconda3/envs/mamba/envs/velocyto2/
gtf=/data/home/zhou3lab/suqianwei/project/seahorse/13.singcell/Ha_testis/Ha_anno_process/updateGSDF/Ha.f2.gtf

velocyto run10x Ha_testis_v9anno $gtf -@4
