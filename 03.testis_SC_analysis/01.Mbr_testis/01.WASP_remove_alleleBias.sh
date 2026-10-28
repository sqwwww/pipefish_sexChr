#!/bin/bash
#SBATCH --job-name=allele
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=5G
#SBATCH --output=%j.out
#SBATCH --error=%j.err


#source /share/home/zhou3lab/suqianwei/miniconda3/bin/activate /share/home/zhou3lab/suqianwei/miniconda3/envs/star2.7.1a
#STAR \
#--genomeDir /data/home/zhou3lab/suqianwei/project/seahorse/13.singcell/Mbr_testis/ref2 \
#--readFilesIn \
#    /data/home/zhou3lab/suqianwei/project/seahorse/01.data/10.singleCell/Mbr_testis/0402_1other/combined_data/cDNA/Mbtcombined_S1_L001_1.fq.gz \
#    /data/home/zhou3lab/suqianwei/project/seahorse/01.data/10.singleCell/Mbr_testis/0402_1other/combined_data/cDNA/Mbtcombined_S1_L001_2.fq.gz \
#--readFilesCommand zcat \
#--runThreadN 40 \
#--outFileNamePrefix ./Mbr_testis_STARmap/map1 \
#--outSAMtype BAM SortedByCoordinate \
#--limitOutSJcollapsed 10000000 \
#--limitIObufferSize 1500000000 



#source /share/home/zhou3lab/suqianwei/miniconda3/bin/activate /share/home/zhou3lab/suqianwei/miniconda3/envs/WASP
#export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
#module load samtools-1.11
## get SNPs from VCF files:
##/data/home/zhou3lab/suqianwei/tools2/WASP-master/mapping/extract_vcf_snps.sh inputVCF processedVCF
#python /data/home/zhou3lab/suqianwei/tools2/WASP-master/mapping/find_intersecting_snps.py \
#       --is_paired_end \
#       --is_sorted \
#       --output_dir outBAM \
#       --snp_dir processedVCF \
#       inputBAM/map1.bam

#source /share/home/zhou3lab/suqianwei/miniconda3/bin/activate /share/home/zhou3lab/suqianwei/miniconda3/envs/star2.7.1a
#STAR \
#--genomeDir /data/home/zhou3lab/suqianwei/project/seahorse/13.singcell/Mbr_testis/ref2 \
#--readFilesIn \
#    outBAM/map1.remap.fq1.gz \
#    outBAM/map1.remap.fq2.gz \
#--readFilesCommand zcat \
#--runThreadN 4 \
#--outFileNamePrefix map2/map2 \
#--outSAMtype BAM SortedByCoordinate \
#--limitOutSJcollapsed 10000000 \
#--limitIObufferSize 1500000000



source /share/home/zhou3lab/suqianwei/miniconda3/bin/activate /share/home/zhou3lab/suqianwei/miniconda3/envs/WASP
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
module load samtools-1.11
python /data/home/zhou3lab/suqianwei/tools2/WASP-master/mapping/filter_remapped_reads.py \
  outBAM/map1.to.remap.bam \
  map2/map2Aligned.sortedByCoord.out.bam \
  filter_remapped_reads/keep_afterFilt.bam

#samtools sort -@ 2 -o filter_remapped_reads/${id}.keep_afterFilt.sort.bam filter_remapped_reads/${id}.keep_afterFilt.bam
#python /data/home/zhou3lab/suqianwei/tools2/WASP-master/mapping/rmdup_pe.py filter_remapped_reads/${id}.keep_afterFilt.sort.bam filter_remapped_reads/${id}.keep_afterFilt.rmdup.bam

#export PATH=/data/home/zhou3lab/suqianwei/tools/gatk-4.3.0.0:$PATH
#fa=../HJ.genome.masked.chr.fa
#gatk ASEReadCounter -R $fa -I filter_remapped_reads/${id}.keep_afterFilt.sort.bam -V ../all_heterM_homoF_snp.m3.vcf -O ${id}.AS.v3.table --minimum-mapping-quality 20 --min-base-quality 5 --min-depth-of-non-filtered-base 10 --read-filter MappingQualityReadFilter

