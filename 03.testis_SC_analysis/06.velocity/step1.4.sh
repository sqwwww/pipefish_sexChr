#!/bin/bash
#SBATCH -J velocyto
#SBATCH -p cpu
#SBATCH -N 1
#SBATCH --mem=30G
#SBATCH --ntasks-per-node=6
#SBATCH --output=%J.out
#SBATCH --error=%J.err




source /share/home/zhou3lab/suqianwei/miniconda3/bin/activate ~/miniconda3/envs/mamba/envs/velocyto2/
gtf=/data/home/zhou3lab/suqianwei/project/seahorse/13.singcell/Mbr_testis/gtfProcess/HJ.addMito.gtf


velocyto run -b output/filter_matrix/barcodes.tsv.gz -o output out.velocyto.bam $gtf -@4 --samtools-memory 4000
