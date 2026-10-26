#!/bin/bash
#SBATCH --job-name=polish
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=25
#SBATCH --mem=180G
#SBATCH --output=polished1.o
#SBATCH --error=polished1.e

/data/home/zhou3lab/suqianwei/tools2/meryl-1.3/bin/meryl count k=15 output merylDB HF_raw.fasta
/data/home/zhou3lab/suqianwei/tools2/meryl-1.3/bin/meryl print greater-than distinct=0.9998 merylDB > repetitive_k15.txt
/share/home/zhou3lab/zhangxinpei/bioapp/miniconda3/bin/winnowmap -t 25 -W repetitive_k15.txt -ax map-pb HF_raw.fasta HF.filt.fastq.gz | samtools view -@ 25 -bS - > hifi.map.bam
samtools sort -@ 25 -o hifi.map.sort.bam hifi.map.bam
samtools index -@ 25 hifi.map.sort.bam
rm hifi.map.bam
/data/home/zhou3lab/wangyifeng/wyf/tools/yak/yak/yak count -o k31.yak -k 31 -b 37 -t 18 <(cat HF_female_DNA_*.fq) <(zcat HF_male_DNA_*.fq.gz)
/data/home/zhou3lab/wangyifeng/wyf/tools/yak/yak/yak count -o k21.yak -k 21 -b 37 -t 18 <(cat HF_female_DNA_*.fq) <(zcat HF_male_DNA_*.fq.gz)


/share/home/zhou3lab/zhangxinpei/bioapp/miniconda3/bin/nextPolish2 -t 25 hifi.map.sort.bam HF_raw.fasta k21.yak k31.yak > HF.polished.round1.fa

/data/home/zhou3lab/wangyifeng/wyf/tools/meryl/meryl-1.3/bin/meryl count k=15 output merylDB HF.polished.round1.fa
/data/home/zhou3lab/wangyifeng/wyf/tools/meryl/meryl-1.3/bin/meryl print greater-than distinct=0.9998 merylDB > repetitive_k15.1.txt
/share/home/zhou3lab/zhangxinpei/bioapp/miniconda3/bin/winnowmap -t 25 -W repetitive_k15.1.txt -ax map-pb HF.polished.round1.fa HF.filt.fastq.gz | samtools view -@ 25 -bS - > hifi.map.1.bam
samtools sort -@ 25 -o hifi.map.sort.1.bam hifi.map.1.bam
samtools index -@ 25 hifi.map.sort.1.bam
#rm hifi.map.bam
/share/home/zhou3lab/zhangxinpei/bioapp/miniconda3/bin/nextPolish2 -t 25 hifi.map.sort.1.bam HF.polished.round1.fa k21.yak k31.yak > HF.polished.round2.fa

/data/home/zhou3lab/wangyifeng/wyf/tools/meryl/meryl-1.3/bin/meryl count k=15 output merylDB HF.polished.round2.fa
/data/home/zhou3lab/wangyifeng/wyf/tools/meryl/meryl-1.3/bin/meryl print greater-than distinct=0.9998 merylDB > repetitive_k15.2.txt
/share/home/zhou3lab/zhangxinpei/bioapp/miniconda3/bin/winnowmap -t 25 -W repetitive_k15.2.txt -ax map-pb HF.polished.round2.fa HF.filt.fastq.gz | samtools view -@ 25 -bS - > hifi.map.2.bam
samtools sort -@ 25 -o hifi.map.sort.2.bam hifi.map.2.bam
samtools index -@ 25 hifi.map.sort.2.bam
#rm hifi.map.bam
/share/home/zhou3lab/zhangxinpei/bioapp/miniconda3/bin/nextPolish2 -t 25 hifi.map.sort.2.bam HF.polished.round2.fa k21.yak k31.yak > HF.polished.round3.fa

