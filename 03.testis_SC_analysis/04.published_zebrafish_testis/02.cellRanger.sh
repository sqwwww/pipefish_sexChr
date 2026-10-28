#!/bin/bash
#SBATCH -J sc
#SBATCH -p short
#SBATCH -N 1
#SBATCH --mem=100G
#SBATCH --ntasks-per-node=30
#SBATCH --output=%J.out
#SBATCH --error=%J.err


cellranger=/data/home/zhou3lab/caicheng/software/cellranger-7.2.0/cellranger
ref_fasta=/data/home/zhou3lab/suqianwei/project/seahorse/02.ref/21.otherSyng/zipfile/zebrafish/Lawson/danRer11.primary.fa
gtf=/data/home/zhou3lab/suqianwei/project/seahorse/02.ref/21.otherSyng/zipfile/zebrafish/Lawson/v4.3.2.gtf
fq_path=../01.data/02.rename/
date
#$cellranger mkref --genome=zbf --fasta=$ref_fasta --genes=$gtf --nthreads=24
date
$cellranger count --id=zbfTestis --transcriptome=zbf --fastqs=$fq_path --sample=zbfTestis
date






