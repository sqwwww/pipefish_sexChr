#!/bin/bash
#SBATCH -J f1
#SBATCH -p cpu
#SBATCH -N 1
#SBATCH --mem=10G
#SBATCH --ntasks-per-node=2
#SBATCH --output=%J.out
#SBATCH --error=%J.err


removeList=../../filter_remapped_reads/remove_readsName.txt

sbatch filter_single_fastq.sh /data/home/zhou3lab/suqianwei/project/seahorse/01.data/10.singleCell/Ha_testis/02.rename/Ha_testis_S1_L001_R1_001.fastq.gz Ha_testisFilt_S1_L001_R1_001.fastq.gz $removeList
sbatch filter_single_fastq.sh /data/home/zhou3lab/suqianwei/project/seahorse/01.data/10.singleCell/Ha_testis/02.rename/Ha_testis_S1_L001_R2_001.fastq.gz Ha_testisFilt_S1_L001_R2_001.fastq.gz $removeList

sbatch filter_single_fastq.sh /data/home/zhou3lab/suqianwei/project/seahorse/01.data/10.singleCell/Ha_testis/02.rename/Ha_testis_S1_L002_R1_001.fastq.gz Ha_testisFilt_S1_L002_R1_001.fastq.gz $removeList
sbatch filter_single_fastq.sh /data/home/zhou3lab/suqianwei/project/seahorse/01.data/10.singleCell/Ha_testis/02.rename/Ha_testis_S1_L002_R2_001.fastq.gz Ha_testisFilt_S1_L002_R2_001.fastq.gz $removeList

sbatch filter_single_fastq.sh /data/home/zhou3lab/suqianwei/project/seahorse/01.data/10.singleCell/Ha_testis/02.rename/Ha_testis_S1_L003_R1_001.fastq.gz Ha_testisFilt_S1_L003_R1_001.fastq.gz $removeList
sbatch filter_single_fastq.sh /data/home/zhou3lab/suqianwei/project/seahorse/01.data/10.singleCell/Ha_testis/02.rename/Ha_testis_S1_L003_R2_001.fastq.gz Ha_testisFilt_S1_L003_R2_001.fastq.gz $removeList

sbatch filter_single_fastq.sh /data/home/zhou3lab/suqianwei/project/seahorse/01.data/10.singleCell/Ha_testis/02.rename/Ha_testis_S1_L004_R1_001.fastq.gz Ha_testisFilt_S1_L004_R1_001.fastq.gz $removeList
sbatch filter_single_fastq.sh /data/home/zhou3lab/suqianwei/project/seahorse/01.data/10.singleCell/Ha_testis/02.rename/Ha_testis_S1_L004_R2_001.fastq.gz Ha_testisFilt_S1_L004_R2_001.fastq.gz $removeList
