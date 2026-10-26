#!/bin/bash
#SBATCH --job-name=tgs
#SBATCH --time=300:00:00 # hh:mm:ss
#SBATCH --partition=cpu
#SBATCH -N 1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=300G
#SBATCH --output=%j.out
#SBATCH --error=%j.err



source /share/home/zhou3lab/suqianwei/miniconda3/bin/activate tgsgapcloser

tgsgapcloser  \
        --scaff  ../HF.ragtag.fasta \
        --reads  /data/home/zhou3lab/suqianwei/project/seahorse/01.data/07.hifi/02.clean/HF.filt.fastq.gz \
        --output HF_nocorrection \
	--ngs /data/home/zhou3lab/suqianwei/project/seahorse/01.data/01.T7clean/02.hf/numName/hf*fq.gz \
        --tgstype pb --thread 4 \
	--pilon /share/home/zhou3lab/suqianwei/miniconda3/envs/tgsgapcloser/bin/pilon \
	--racon /share/home/zhou3lab/suqianwei/miniconda3/envs/tgsgapcloser/bin/racon \
	--samtools /share/home/zhou3lab/suqianwei/miniconda3/envs/tgsgapcloser/bin/samtools\
	--java /share/home/zhou3lab/suqianwei/miniconda3/envs/tgsgapcloser/bin/java \
	--ne
