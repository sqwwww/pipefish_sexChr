#!/bin/bash
#SBATCH --job-name=PB.hifiasm
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=100000
#SBATCH --output=PB.hifiasm.o
#SBATCH --error=PB.hifiasm.e

hic_fq1=/data/home/zhou3lab/suqianwei/project/seahorse/01.data/08.hic/02.HJ/HJ-M_R1.fq.gz
hic_fq2=/data/home/zhou3lab/suqianwei/project/seahorse/01.data/08.hic/02.HJ/HJ-M_R2.fq.gz
hifiReads=/data/home/zhou3lab/suqianwei/project/seahorse/01.data/07.hifi/02.clean/HJ.filt.fastq.gz
/data/home/zhou3lab/suqianwei/tools2/hifiasm-master/hifiasm -o HJ.hifiasm -t 16 --h1 ${hic_fq1} --h2 ${hic_fq2} ${hifiReads}
awk '/^S/{print ">"$2;print $3}' HJ.hifiasm.hic.p_ctg.gfa > HJ.p_ctg.fa
