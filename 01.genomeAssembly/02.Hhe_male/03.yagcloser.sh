#!/bin/bash

#SBATCH --job-name=yagcloser
#SBATCH --partition=short
#SBATCH -N 1
#SBATCH --ntasks-per-node=12
#SBATCH --mem=50g
#SBATCH --output=yagcloser.o
#SBATCH --error=yagcloser.e
#conda activate yagcloser

source /share/home/zhou3lab/suqianwei/miniconda3/bin/activate yagcloser
module load samtools-1.11

reads=/data/home/zhou3lab/suqianwei/project/seahorse/01.data/07.hifi/02.clean/HF.filt.fastq.gz
fa=../03.ragtag/Cin_outgroup_hifi/ragtag.scaffold.fasta

#/data/home/zhou3lab/suqianwei/tools2/minimap2-2.26_x64-linux/minimap2 -t 12 --secondary=yes -ax map-hifi $fa $reads | samtools view -@ 12 -hSb -q 20 - | samtools sort -@ 12 -m 4G - > aln.s.bam
#samtools index -@ 8 aln.s.bam
#/data/home/zhou3lab/suqianwei/tools2/asset/bin/detgaps $fa > gaps.bed
python /data/home/zhou3lab/suqianwei/tools2/yagcloser-master/yagcloser.py -g $fa -a aln.s.bam -b gaps.bed -o yagcloser_output -f 40 -mins 2 -s HF
#python /data/home/zhou3lab/suqianwei/tools2/yagcloser-master/scripts/update_assembly_edits_and_breaks.py -i $fa -o final.eidt.fasta -e yagcloser_output/HF.edits.txt
