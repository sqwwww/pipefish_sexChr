#!/bin/bash
#SBATCH --job-name=tgs
#SBATCH --time=300:00:00 # hh:mm:ss
#SBATCH --partition=cpu
#SBATCH -N 1
#SBATCH --ntasks-per-node=12
#SBATCH --mem=300G
#SBATCH --output=%j.out
#SBATCH --error=%j.err



source /share/home/zhou3lab/suqianwei/miniconda3/bin/activate tgsgapcloser

tgsgapcloser  \
        --scaff  /data/home/zhou3lab/suqianwei/project/seahorse/09.genomeAsemb/HJ/06.gapfill/01.yagcloser/final.eidt.fasta \
        --reads  /data/home/zhou3lab/suqianwei/project/seahorse/01.data/07.hifi/02.clean/HJ.filt.fastq.gz \
        --output HJ_nocorrection \
        --tgstype pb --thread 12 \
	--samtools /share/home/zhou3lab/suqianwei/miniconda3/envs/tgsgapcloser/bin/samtools\
	--java /share/home/zhou3lab/suqianwei/miniconda3/envs/tgsgapcloser/bin/java \
	--ne
