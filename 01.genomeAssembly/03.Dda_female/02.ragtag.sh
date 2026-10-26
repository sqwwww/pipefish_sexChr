#!/bin/bash
#SBATCH --job-name=ragtag
#SBATCH --partition=short
#SBATCH -N 1
#SBATCH --mem=30G
#SBATCH --ntasks-per-node=12
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH --time=12:00:00


source /share/home/zhou3lab/suqianwei/miniconda3/bin/activate ragtag
ref=/data/home/zhou3lab/suqianwei/project/seahorse/02.ref/21.otherSyng/genome_fasta/Dunckerocampus_dactyliophorus.genome.fasta
qry=../01.hifiasm/BJ.p_ctg.fasta
ragtag.py scaffold $ref $qry -t 12

source /share/home/zhou3lab/suqianwei/miniconda3/bin/activate quast
quast.py -t 10 ragtag_output/ragtag.scaffold.fasta

