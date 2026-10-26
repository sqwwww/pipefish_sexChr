#!/bin/bash
#SBATCH --job-name=HF.hifiasm
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=20G
#SBATCH --output=HF.hifiasm.o
#SBATCH --error=HF.hifiasm.e

hifiReads=/data/home/zhou3lab/suqianwei/project/seahorse/01.data/07.hifi/02.clean/HF.filt.fastq.gz
/data/home/zhou3lab/suqianwei/tools2/hifiasm-master/hifiasm -o HF.hifiasm -t 16 ${hifiReads}
awk '/^S/{print ">"$2;print $3}' HF.hifiasm.bp.p_ctg.gfa > HF.p_ctg.fa
