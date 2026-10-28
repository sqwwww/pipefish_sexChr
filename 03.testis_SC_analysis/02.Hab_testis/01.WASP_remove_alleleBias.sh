#!/bin/bash
#SBATCH --job-name=allele
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=5G
#SBATCH --output=%j.out
#SBATCH --error=%j.err


#source /share/home/zhou3lab/suqianwei/miniconda3/bin/activate /share/home/zhou3lab/suqianwei/miniconda3/envs/star2.7.1a
#
#STAR \
#  --genomeDir ../HaDB_v6/star/ \
#  --readFilesIn inputBAM/Ha_testis_combined_R2.fastq.gz \
#  --readFilesCommand zcat \
#  --runThreadN 40 \
#  --outFileNamePrefix inputBAM/map1 \
#  --outFilterScoreMin 30 \
#  --outFilterMultimapNmax 1 \
#  --alignIntronMax 1000000 \
#  --outSAMtype BAM SortedByCoordinate \
#  --outSAMattributes NH HI AS nM
#




#source /share/home/zhou3lab/suqianwei/miniconda3/bin/activate /share/home/zhou3lab/suqianwei/miniconda3/envs/WASP
#export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
#module load samtools-1.11
### get SNPs from VCF files:
##/data/home/zhou3lab/suqianwei/tools2/WASP-master/mapping/extract_vcf_snps.sh inputVCF processedVCF
##samtools index inputBAM/${id}.rg.bam
#python /data/home/zhou3lab/suqianwei/tools2/WASP-master/mapping/find_intersecting_snps.py \
#       --is_sorted \
#       --output_dir outBAM \
#       --snp_dir processedVCF \
#       inputBAM/map1Aligned.sortedByCoord.out.bam


## remap selected processed reads
#source /share/home/zhou3lab/suqianwei/miniconda3/bin/activate /share/home/zhou3lab/suqianwei/miniconda3/envs/star2.7.1a
#
#STAR \
#  --genomeDir ../HaDB_v6/star/ \
#  --readFilesIn outBAM/map1Aligned.sortedByCoord.out.remap.fq.gz \
#  --readFilesCommand zcat \
#  --runThreadN 40 \
#  --outFileNamePrefix map2/map2 \
#  --outFilterScoreMin 30 \
#  --outFilterMultimapNmax 1 \
#  --alignIntronMax 1000000 \
#  --outSAMtype BAM SortedByCoordinate \
#  --outSAMattributes NH HI AS nM




source /share/home/zhou3lab/suqianwei/miniconda3/bin/activate /share/home/zhou3lab/suqianwei/miniconda3/envs/WASP
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
module load samtools-1.11
python /data/home/zhou3lab/suqianwei/tools2/WASP-master/mapping/filter_remapped_reads.py \
  outBAM/map1Aligned.sortedByCoord.out.to.remap.bam \
  map2/map2Aligned.sortedByCoord.out.bam \
  filter_remapped_reads/keep_afterFilt.bam

#samtools sort -@ 2 -o filter_remapped_reads/${id}.keep_afterFilt.sort.bam filter_remapped_reads/${id}.keep_afterFilt.bam
#python /data/home/zhou3lab/suqianwei/tools2/WASP-master/mapping/rmdup_pe.py filter_remapped_reads/${id}.keep_afterFilt.sort.bam filter_remapped_reads/${id}.keep_afterFilt.rmdup.bam
#
#export PATH=/data/home/zhou3lab/suqianwei/tools/gatk-4.3.0.0:$PATH
#fa=${ref_fa}
##gatk CreateSequenceDictionary -R $fa -O ${fa%fasta}dict
##sed 's/Chr//' ../all_heterM_homoF_snp.m4.vcf > all_heterM_homoF_snp.m5.vcf
##gatk IndexFeatureFile -I all_heterM_homoF_snp.m5.vcf
#samtools view -h filter_remapped_reads/${id}.keep_afterFilt.sort.bam | sed 's/Chr//g' > filter_remapped_reads/${id}.keep_afterFilt.v2.sort.bam
#gatk ASEReadCounter -R $fa -I filter_remapped_reads/${id}.keep_afterFilt.v2.sort.bam -V all_heterM_homoF_snp.m5.vcf -O ${id}.AS.v3.table --minimum-mapping-quality 20 --min-base-quality 5 --min-depth-of-non-filtered-base 0 --read-filter MappingQualityReadFilter
#
