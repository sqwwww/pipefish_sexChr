#!/bin/bash
#SBATCH --job-name=ragtag
#SBATCH --partition=short
#SBATCH -N 1
#SBATCH --mem=100G
#SBATCH --ntasks-per-node=16
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH --time=12:00:00

## prepare
#cp /data/home/zhou3lab/suqianwei/project/seahorse/09.genomeAsemb/20.stflr/hf2/hf222_output.fasta.gz hf_stlfr.fasta.gz
#gunzip hf_stlfr.fasta.gz
#ln -s /data/home/zhou3lab/suqianwei/project/seahorse/02.ref/03.seahorse_Ha/Ha_genome.chromosome.fasta ./

source /share/home/zhou3lab/suqianwei/miniconda3/bin/activate ragtag
ref=/data/home/zhou3lab/suqianwei/project/seahorse/02.ref/21.otherSyng/genome_fasta/Corythoichthys_intestinalis.genome.fasta
qry=../01.hifiasm/HF.p_ctg.fasta
ragtag.py scaffold $ref $qry -t 16 -o Cin_outgroup_hifi

#source /share/home/zhou3lab/suqianwei/miniconda3/bin/activate quast
#quast.py -t 10 Ddac_outgroup/ragtag.scaffold.fasta

