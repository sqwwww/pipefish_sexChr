#!/bin/bash
#SBATCH --job-name=merqury
#SBATCH --time=12:00:00 # hh:mm:ss
#SBATCH --partition=short
#SBATCH -N 1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=70G
#SBATCH --output=%j.out
#SBATCH --error=%j.err


#k=19
## 1. Build meryl dbs on each input, could be run on different nodes
#/data/home/zhou3lab/suqianwei/tools2/meryl-1.3/bin/meryl k=$k count output HF_female_DNA_1.meryl HF_female_DNA_1.fq.gz
#/data/home/zhou3lab/suqianwei/tools2/meryl-1.3/bin/meryl k=$k count output HF_female_DNA_2.meryl HF_female_DNA_2.fq.gz
#
#/data/home/zhou3lab/suqianwei/tools2/meryl-1.3/bin/meryl k=$k count output HF_male_DNA_1.meryl HF_male_DNA_1.fq.gz
#/data/home/zhou3lab/suqianwei/tools2/meryl-1.3/bin/meryl k=$k count output HF_male_DNA_2.meryl HF_male_DNA_2.fq.gz
## 2. Merge
#/data/home/zhou3lab/suqianwei/tools2/meryl-1.3/bin/meryl union-sum output combined.meryl HF_*.meryl
#

source /share/home/zhou3lab/suqianwei/miniconda3/bin/activate /share/home/zhou3lab/suqianwei/miniconda3/envs/R4
export MERQURY=/data/home/zhou3lab/suqianwei/tools2/merqury-master
export PATH=/data/home/zhou3lab/suqianwei/tools2/meryl-1.4.1/bin/:$PATH
# 3. merqury
/data/home/zhou3lab/suqianwei/tools2/merqury-master/merqury.sh combined.meryl HF.polished.round3.fa merqury_test1

