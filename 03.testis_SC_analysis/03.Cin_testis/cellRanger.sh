#!/bin/bash
#SBATCH -J sc
#SBATCH -p cpu
#SBATCH -N 1
#SBATCH --mem=100G
#SBATCH --ntasks-per-node=24
#SBATCH --output=%J.out
#SBATCH --error=%J.err


cellranger=/data/home/zhou3lab/caicheng/software/cellranger-7.2.0/cellranger
ref_fasta=ref_file/Cint.combined_genome_mito.fa
gtf=ref_file/Cint.combined_genome_mito.gtf
fq_path=/data/home/zhou3lab/suqianwei/project/seahorse/01.data/10.singleCell/Cint_testis
#date
#$cellranger mkref --genome=CintDB_v1 --fasta=$ref_fasta --genes=$gtf --nthreads=24
date
# v1
#$cellranger count --id=CintTestis_v1 --transcriptome=CintDB_v1 --fastqs=$fq_path --sample=PA_mt-3,PA_mt-2,PA_mt-1 --localcores 24 

# v2
#$cellranger count --id=CintTestis_v2 --transcriptome=CintDB_v1 --fastqs=$fq_path --sample=PA_mt-3,PA_mt-2,PA_mt-1 --localcores 24 --force-cells=4000

# v3
$cellranger count --id=CintTestis_v3 --transcriptome=CintDB_v1 --fastqs=$fq_path --sample=PA_mt-3,PA_mt-2,PA_mt-1 --localcores 24 --force-cells=3000

date





