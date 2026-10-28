#!/bin/bash
#SBATCH -J sc
#SBATCH -p cpu
#SBATCH -N 1
#SBATCH --mem=30G
#SBATCH --ntasks-per-node=24
#SBATCH --output=%J.out
#SBATCH --error=%J.err


cellranger=/data/home/zhou3lab/caicheng/software/cellranger-7.2.0/cellranger
ref_fasta=/data/home/zhou3lab/suqianwei/project/seahorse/13.singcell/Ha_testis/Ha_ChrGenome_Mito.fasta
gtf=/data/home/zhou3lab/suqianwei/project/seahorse/13.singcell/Ha_testis/Ha_anno_process/updateGSDF/Ha.f2.gtf
fq_path=/data/home/zhou3lab/suqianwei/project/seahorse/13.singcell/Ha_testis/removeAlleleBias/cellRanger/01.processedFq/
#date
#$cellranger mkref --genome=HaDB_v6 --fasta=$ref_fasta --genes=$gtf --nthreads=8
date
$cellranger count --id=Ha_testis_v9anno --transcriptome=/data/home/zhou3lab/suqianwei/project/seahorse/13.singcell/Ha_testis/HaDB_v6 --fastqs=$fq_path --sample=Ha_testisFilt --localcores 24
date





