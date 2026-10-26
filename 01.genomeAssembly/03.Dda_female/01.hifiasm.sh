#!/bin/bash
#SBATCH --job-name=BJ.hifiasm
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=20
#SBATCH --mem=20G
#SBATCH --output=BJ.hifiasm.o
#SBATCH --error=BJ.hifiasm.e

hifiReads=/data/home/zhou3lab/suqianwei/project/seahorse/01.data/07.hifi/02.clean/BJ.filt.fastq.gz
/data/home/zhou3lab/suqianwei/tools2/hifiasm-master/hifiasm -o HJ.hifiasm -t 20 ${hifiReads}
