#!/bin/bash
#SBATCH -J sc
#SBATCH -p cpu
#SBATCH -N 1
#SBATCH --mem=10G
#SBATCH --ntasks-per-node=4
#SBATCH --output=%J.out
#SBATCH --error=%J.err



# after removeing reads that would cause allele-mapping bias, I rerun the dnbc4tools
dnbc4tools=/data/home/zhou3lab/suqianwei/tools2/dnbc4tools2.1.3/dnbc4tools

gtf=gtfProcess/HJ.addMito.gtf
fa=gtfProcess/HJ.addMito.fa
#$dnbc4tools tools mkgtf --action check --ingtf $gtf --output corrected.gtf
#$dnbc4tools rna mkref --ingtf corrected.gtf --fasta $fa --threads 20 --species Mbr --genomeDir ref
#$dnbc4tools rna mkref --ingtf gtfProcess/geneExt/HJ/result.gtf --fasta $fa --threads 20 --species Mbr --genomeDir ref2

#$dnbc4tools rna run \
#		--name Mbr_testis3 \
#		--cDNAfastq1 /data/home/zhou3lab/suqianwei/project/seahorse/13.singcell/Mbr_testis/removeAlleleBias/cellRanger/01.processedFq/cDNA/MbtcombinedFilt_S1_L001_1.fq.gz \
#		--cDNAfastq2 /data/home/zhou3lab/suqianwei/project/seahorse/13.singcell/Mbr_testis/removeAlleleBias/cellRanger/01.processedFq/cDNA/MbtcombinedFilt_S1_L001_2.fq.gz \
#		--oligofastq1 /data/home/zhou3lab/suqianwei/project/seahorse/01.data/10.singleCell/Mbr_testis/0402_1other/combined_data/oligo/Mbtcombined_S1_L001_1.fq.gz \
#		--oligofastq2 /data/home/zhou3lab/suqianwei/project/seahorse/01.data/10.singleCell/Mbr_testis/0402_1other/combined_data/oligo/Mbtcombined_S1_L001_2.fq.gz \
#		--genomeDir /data/home/zhou3lab/suqianwei/project/seahorse/13.singcell/Mbr_testis/ref2 \
#		--threads 20

$dnbc4tools tools changetag --inbam Mbr_testis3/output/anno_decon_sorted.bam --outbam out.velocyto.bam

