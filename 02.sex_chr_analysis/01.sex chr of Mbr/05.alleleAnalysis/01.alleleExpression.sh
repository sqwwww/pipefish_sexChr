#!/bin/bash
#SBATCH --job-name=allele
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=20G
#SBATCH --output=%j.out
#SBATCH --error=%j.err


source /share/home/zhou3lab/suqianwei/miniconda3/bin/activate /share/home/zhou3lab/suqianwei/miniconda3/envs/WASP
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
module load samtools-1.11
### get SNPs from VCF files:
/data/home/zhou3lab/suqianwei/tools2/WASP-master/mapping/extract_vcf_snps.sh inputVCF processedVCF
id=Mba_Mp
samtools index inputBAM/${id}.rg.bam
python /data/home/zhou3lab/suqianwei/tools2/WASP-master/mapping/find_intersecting_snps.py \
       --is_paired_end \
       --is_sorted \
       --output_dir outBAM \
       --snp_dir processedVCF \
       inputBAM/${id}.rg.bam


source /share/home/zhou3lab/suqianwei/miniconda3/bin/activate /share/home/zhou3lab/suqianwei/miniconda3/envs/transcriptome
ref_fa=/data/home/zhou3lab/suqianwei/project/seahorse/09.genomeAsemb/HJ/08.removeBlankRegion/final2.chr.fasta
hisat2-build -p 20 $ref_fa hisat2Index/HJ
gtf=/data/home/zhou3lab/suqianwei/project/seahorse/06.genomeAnno/06.braker/08.HJ_hifi_final/braker/braker.gtf
fq1=outBAM/${id}.rg.remap.fq1.gz
fq2=outBAM/${id}.rg.remap.fq2.gz
ref=../../01.bam/01.HJ/hisat2Index/HJ
hisat2 -x $ref -1 $fq1 -2 $fq2 --phred33 -p 40 | samtools view -q 30 -b -@ 10 | samtools sort -@ 6 -m 1g -o map2/${id}.remap.sorted.bam
samtools index map2/${id}.remap.sorted.bam

source /share/home/zhou3lab/suqianwei/miniconda3/bin/activate /share/home/zhou3lab/suqianwei/miniconda3/envs/WASP
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
module load samtools-1.11
python /data/home/zhou3lab/suqianwei/tools2/WASP-master/mapping/filter_remapped_reads.py \
  outBAM/${id}.rg.to.remap.bam \
  map2/${id}.remap.sorted.bam \
  filter_remapped_reads/${id}.keep_afterFilt.bam

samtools sort -@ 2 -o filter_remapped_reads/${id}.keep_afterFilt.sort.bam filter_remapped_reads/${id}.keep_afterFilt.bam
python /data/home/zhou3lab/suqianwei/tools2/WASP-master/mapping/rmdup_pe.py filter_remapped_reads/${id}.keep_afterFilt.sort.bam filter_remapped_reads/${id}.keep_afterFilt.rmdup.bam

export PATH=/data/home/zhou3lab/suqianwei/tools/gatk-4.3.0.0:$PATH
fa=../HJ.genome.masked.chr.fa
gatk ASEReadCounter -R $fa -I filter_remapped_reads/${id}.keep_afterFilt.sort.bam -V ../all_heterM_homoF_snp.m3.vcf -O ${id}.AS.v3.table --minimum-mapping-quality 20 --min-base-quality 5 --min-depth-of-non-filtered-base 10 --read-filter MappingQualityReadFilter

