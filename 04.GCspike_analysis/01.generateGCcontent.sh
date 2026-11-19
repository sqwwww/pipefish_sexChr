#!/bin/bash
#SBATCH --job-name=repeatdensity
#SBATCH --time=12:00:00 # hh:mm:ss
#SBATCH --partition=mem
#SBATCH -N 1
#SBATCH --ntasks-per-node=2
#SBATCH --output=%j.out
#SBATCH --error=%j.err


module load bedtools-2.29.2


bedtools nuc -fi ../../Cin.numChr.fasta -bed Cin_chrs.100k | awk -v OFS='\t' '{print $1,$2,$3,($7+$8)/($12-$11-$10)}' | sed '1,1d' > Cin.GCcontent.100k.bed
bedtools nuc -fi /share/home/zhou3lab/suqianwei/project/seahorse/02.ref/01.ringedPipe/ringedPipefish.chromosome.fasta -bed BJM_chrs.100k | awk -v OFS='\t' '{print $1,$2,$3,($7+$8)/($12-$11-$10)}' | sed '1,1d' > BJ.GCcontent.100k.bed
bedtools nuc -fi ../../HF.genome.chr.fasta -bed HF_chrs.100k | awk -v OFS='\t' '{print $1,$2,$3,($7+$8)/($12-$11-$10)}' | sed '1,1d' > HF.GCcontent.100k.bed
bedtools nuc -fi ../../HJ.genome.chr.fasta -bed HJ_chrs.100k | awk -v OFS='\t' '{print $1,$2,$3,($7+$8)/($12-$11-$10)}' | sed '1,1d' > HJ.GCcontent.100k.bed
# 5th column was GC content of the window
