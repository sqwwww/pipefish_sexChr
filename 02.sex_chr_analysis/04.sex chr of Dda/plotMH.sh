#!/bin/bash
#SBATCH --job-name=Manhattan
#SBATCH --partition=short
#SBATCH -N 1
#SBATCH --mem=20g
#SBATCH --ntasks-per-node=1
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH --time=6:00:00

source /share/home/zhou3lab/suqianwei/miniconda3/bin/activate /share/home/zhou3lab/suqianwei/miniconda3/envs/R4
#Rscript plotMH.R dataForManhattan.aIBS.txt BJM_hifi_aIBS
Rscript plotMH.R dataForManhattan.aBN.txt BJM_hifi_aBN2
