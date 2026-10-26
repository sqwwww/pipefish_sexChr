#!/bin/bash
#SBATCH --job-name=DNA-map
#SBATCH --time=100:00:00 # hh:mm:ss
#SBATCH --partition=cpu
#SBATCH -N 1
#SBATCH --ntasks-per-node=40
#SBATCH --mem=100G
#SBATCH --output=%j.out
#SBATCH --error=%j.err

module load bowtie2-2.4.1
module load samtools-1.11
thread=40
thread2=9



refFa=BJ.genome.masked.chr.fa
bowtie2-build --threads 10 $refFa bowtie2Index/BJmask
bowtie_db=bowtie2Index/BJmask

# used by haplotype caller
export PATH=/data/home/zhou3lab/suqianwei/tools/gatk-4.3.0.0:$PATH
samtools faidx $refFa
gatk CreateSequenceDictionary \
    -R $refFa \
    -O ${refFa%fa}dict

cat sample.list| while read id; do
fq1=/data/home/zhou3lab/suqianwei/project/seahorse/01.data/01.T7clean/04.bj/numName/bj_${id}_1.clean.fq.gz
fq2=/data/home/zhou3lab/suqianwei/project/seahorse/01.data/01.T7clean/04.bj/numName/bj_${id}_2.clean.fq.gz
echo $id
bowtie2 -p ${thread} -x $bowtie_db -1 $fq1 -2 $fq2 | samtools view -b -q 30 -@ ${thread} | samtools sort -@ ${thread2} -m 9g -o ${id}.bam
samtools index ${id}.bam

done


cat sample.list | while read id; do
echo "#!/bin/bash
#SBATCH --job-name=gatk
#SBATCH --time=300:00:00 # hh:mm:ss
#SBATCH --partition=cpu
#SBATCH -N 1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=110G
#SBATCH --output=${id}.haplotypeCaller.out
#SBATCH --error=${id}.haplotypeCaller.err

module load samtools-1.11

export PATH=/data/home/zhou3lab/suqianwei/tools/gatk-4.3.0.0:\$PATH
gatk MarkDuplicates -I ${id}.bam -O $id.dedup.bam -M $id.m
samtools index $id.dedup.bam
gatk AddOrReplaceReadGroups -I $id.dedup.bam -O $id.dedup.rg.bam --RGID 4$id --RGLB lib$id --RGPL illumina --RGPU unit1 --RGSM $id 

samtools index ${id}.dedup.rg.bam

gatk --java-options "-Xmx100G" HaplotypeCaller -R $refFa -I $id.dedup.rg.bam -ERC GVCF -O $id.dedup.rg.gvcf.gz --native-pair-hmm-threads 4
" > ${id}.gatk.sh
sbatch ${id}.gatk.sh
done

