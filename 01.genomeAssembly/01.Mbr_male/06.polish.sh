#!/bin/bash
#SBATCH --job-name=polish
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=25
#SBATCH --mem=180G
#SBATCH --output=%j.out
#SBATCH --error=%j.err

/data/home/zhou3lab/suqianwei/tools2/meryl-1.3/bin/meryl count k=15 output merylDB HJ_raw.fasta
/data/home/zhou3lab/suqianwei/tools2/meryl-1.3/bin/meryl print greater-than distinct=0.9998 merylDB > repetitive_k15.txt
/share/home/zhou3lab/zhangxinpei/bioapp/miniconda3/bin/winnowmap -t 25 -W repetitive_k15.txt -ax map-pb HJ_raw.fasta HJ.filt.fastq.gz | samtools view -@ 25 -bS - > hifi.map.bam
samtools sort -@ 25 -o hifi.map.sort.bam hifi.map.bam
samtools index -@ 25 hifi.map.sort.bam
#rm hifi.map.bam
#/data/home/zhou3lab/wangyifeng/wyf/tools/yak/yak/yak count -o k31.yak -k 31 -b 37 -t 18 <(zcat HJ_female_DNA_*.fq.gz) <(zcat HJ_male_DNA_*.fq.gz)
#/data/home/zhou3lab/wangyifeng/wyf/tools/yak/yak/yak count -o k21.yak -k 21 -b 37 -t 18 <(zcat HJ_female_DNA_*.fq.gz) <(zcat HJ_male_DNA_*.fq.gz)


/share/home/zhou3lab/zhangxinpei/bioapp/miniconda3/bin/nextPolish2 -t 25 hifi.map.sort.bam HJ_raw.fasta k21.yak k31.yak > HJ.polished.round1.fa

/data/home/zhou3lab/wangyifeng/wyf/tools/meryl/meryl-1.3/bin/meryl count k=15 output merylDB HJ.polished.round1.fa
/data/home/zhou3lab/wangyifeng/wyf/tools/meryl/meryl-1.3/bin/meryl print greater-than distinct=0.9998 merylDB > repetitive_k15.1.txt
/share/home/zhou3lab/zhangxinpei/bioapp/miniconda3/bin/winnowmap -t 25 -W repetitive_k15.1.txt -ax map-pb HJ.polished.round1.fa HJ.filt.fastq.gz | samtools view -@ 25 -bS - > hifi.map.1.bam
samtools sort -@ 25 -o hifi.map.sort.1.bam hifi.map.1.bam
samtools index -@ 25 hifi.map.sort.1.bam
#rm hifi.map.bam
/share/home/zhou3lab/zhangxinpei/bioapp/miniconda3/bin/nextPolish2 -t 25 hifi.map.sort.1.bam HJ.polished.round1.fa k21.yak k31.yak > HJ.polished.round2.fa

/data/home/zhou3lab/wangyifeng/wyf/tools/meryl/meryl-1.3/bin/meryl count k=15 output merylDB HJ.polished.round2.fa
/data/home/zhou3lab/wangyifeng/wyf/tools/meryl/meryl-1.3/bin/meryl print greater-than distinct=0.9998 merylDB > repetitive_k15.2.txt
/share/home/zhou3lab/zhangxinpei/bioapp/miniconda3/bin/winnowmap -t 25 -W repetitive_k15.2.txt -ax map-pb HJ.polished.round2.fa HJ.filt.fastq.gz | samtools view -@ 25 -bS - > hifi.map.2.bam
samtools sort -@ 25 -o hifi.map.sort.2.bam hifi.map.2.bam
samtools index -@ 25 hifi.map.sort.2.bam
#rm hifi.map.bam
/share/home/zhou3lab/zhangxinpei/bioapp/miniconda3/bin/nextPolish2 -t 25 hifi.map.sort.2.bam HJ.polished.round2.fa k21.yak k31.yak > HJ.polished.round3.fa

# sort chrs by length
seqkit sort -r --by-length HJ.polished.round3.fa | seqtk rename - sca_ > HJ.polished.round3.fa.sorted.fa
# change chr name to number
sed -i 's/sca_\([0-9]\) /\1 /' HJ.polished.round3.fa.sorted.fa
sed -i 's/sca_\(1[0-9]\) /\1 /' HJ.polished.round3.fa.sorted.fa 
# the final file called HJ.genome.fasta
ln -s HJ.polished.round3.fa.sorted.fa HJ.genome.fasta
# only chr sequence
seqkit grep -f chrs.txt HJ.genome.fasta > HJ.genome.chr.fasta
