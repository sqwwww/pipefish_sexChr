#!/bin/bash
#SBATCH --job-name=quast
#SBATCH --partition=short
#SBATCH -N 1
#SBATCH --mem=10g
#SBATCH --ntasks-per-node=10
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH --time=12:00:00


source /share/home/zhou3lab/suqianwei/miniconda3/bin/activate quast
quast.py -t 10 BJ_F.genome.chr.fasta

