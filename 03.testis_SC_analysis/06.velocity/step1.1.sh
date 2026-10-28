#!/bin/bash
#SBATCH -J velocity
#SBATCH -p cpu
#SBATCH -N 1
#SBATCH --mem=100G
#SBATCH --ntasks-per-node=15
#SBATCH --output=%J.out
#SBATCH --error=%J.err

source /share/home/zhou3lab/suqianwei/miniconda3/bin/activate  /share/home/zhou3lab/suqianwei/miniconda3/envs/mamba/envs/velocyto2
module load samtools-1.11
velocyto run10x -@10 --samtools-memory 8000 ../Ha_testis/Ha_testis_v9anno ../Ha_testis/Ha_anno_process/updateGSDF/Ha.f2.gtf
