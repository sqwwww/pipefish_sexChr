#!/bin/bash
#SBATCH --job-name=HFcoverage
#SBATCH --partition=short
#SBATCH -N 1
#SBATCH --mem=10G
#SBATCH --ntasks-per-node=4
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH --time=12:00:00


#module load bedtools-2.29.2
#
##01. make window of ref
#bedtools makewindows -g HF.genome.masked.chr.fa.fai -w 10000 > HF.genome.masked.chr.fa.10k.bed
#loc_bed=HF.genome.masked.chr.fa.10k.bed
#repeat_bed=/share/home/zhou3lab/suqianwei/project/seahorse/03.repeat/09.HF_hifi/HF_soft/HF.genome.chr.fasta.out.gff.bed
#
##02. combined male bam and caluelate depth
#samtools merge male.bam M1.bam M2.bam M3.bam M4.bam M5.bam M6.bam M7.bam M8.bam M9.bam
#samtools index male.bam
#samtools depth -a male.bam > male_samtools.depth
#awk '{print $1"\t"$2-1"\t"$2"\t"$3}' male_samtools.depth | sort -k1,1n -k2,2n > male_samtools.depth.bed
#bedtools subtract -a male_samtools.depth.bed -b ${repeat_bed} > male_samtools_subtractRepeat.depth.bed
#bedtools map -a $loc_bed -b male_samtools_subtractRepeat.depth.bed -c 4 -o sum,count,median,mean > male_samtools.10kwindow_Srepeat.depth 
#
##03. for female depth
#samtools merge female.bam F1.bam F2.bam F3.bam F4.bam F5.bam F6.bam F7.bam F8.bam F9.bam
#samtools index female.bam
#samtools depth -a female.bam > female_samtools.depth
#awk '{print $1"\t"$2-1"\t"$2"\t"$3}' female_samtools.depth | sort -k1,1n -k2,2n > female_samtools.depth.bed
#bedtools subtract -a female_samtools.depth.bed -b $repeat_bed > female_samtools_subtractRepeat.depth.bed
#bedtools map -a $loc_bed -b female_samtools_subtractRepeat.depth.bed -c 4 -o sum,count,median,mean > female_samtools.10kwindow_Srepeat.depth
#
# extract chr1 mapping result from male bam
samtools view -b -h male.bam 1 > hifi.mapChr1.sort.bam
