#!/bin/bash
#SBATCH --job-name=busco
#SBATCH --partition=short
#SBATCH -N 1
#SBATCH --mem=50g
#SBATCH --ntasks-per-node=20
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH --time=12:00:00
source /share/home/zhou3lab/suqianwei/miniconda3/bin/activate busco
busco -m genome -i HJ.genome.chr.fasta -o HJ.genome.chr.fasta.busco -l /data/home/zhou3lab/suqianwei/tools2/busco_db/actinopterygii_odb10 -c 20 --offline
