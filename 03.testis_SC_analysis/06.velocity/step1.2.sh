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
velocyto run10x -@10 --samtools-memory 8000 ../Cint_testis/CintTestis_v3/ Cint.combined_genome_mito.m1.gtf
