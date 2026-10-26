#!/bin/bash
#SBATCH --job-name=genotype
#SBATCH --time=12:00:00 # hh:mm:ss
#SBATCH --partition=short
#SBATCH -N 1
#SBATCH --ntasks-per-node=5
#SBATCH --mem=70G
#SBATCH --output=%j.out
#SBATCH --error=%j.err

export PATH=/data/home/zhou3lab/suqianwei/tools/gatk-4.3.0.0:$PATH
ref=HJ.genome.masked.chr.fa
#01.CombineGVCFs
gatk --java-options "-Xmx70G" CombineGVCFs \
-R $ref \
-V F1.dedup.rg.gvcf.gz \
-V F2.dedup.rg.gvcf.gz \
-V F3.dedup.rg.gvcf.gz \
-V F4.dedup.rg.gvcf.gz \
-V F5.dedup.rg.gvcf.gz \
-V F6.dedup.rg.gvcf.gz \
-V F7.dedup.rg.gvcf.gz \
-V F8.dedup.rg.gvcf.gz \
-V F9.dedup.rg.gvcf.gz \
-V F10.dedup.rg.gvcf.gz \
-V M1.dedup.rg.gvcf.gz \
-V M2.dedup.rg.gvcf.gz \
-V M3.dedup.rg.gvcf.gz \
-V M4.dedup.rg.gvcf.gz \
-V M5.dedup.rg.gvcf.gz \
-V M6.dedup.rg.gvcf.gz \
-V M7.dedup.rg.gvcf.gz \
-V M8.dedup.rg.gvcf.gz \
-V M9.dedup.rg.gvcf.gz \
-V M10.dedup.rg.gvcf.gz \
-O combined.raw.gvcf.gz
##02.genotype
gatk --java-options "-Xmx70G" GenotypeGVCFs -R $ref -V combined.raw.gvcf.gz -O All.raw.vcf.gz
#03.filter
gatk --java-options "-Xmx70G" SelectVariants -R $ref -V All.raw.vcf.gz -select-type SNP -O ref.snp.vcf
gatk --java-options "-Xmx70G" VariantFiltration -R $ref -V ref.snp.vcf -filter " QD < 2.0 || FS > 60.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 3.0 || MQ < 40.0 " --filter-name "my_snp_filter" -O ref.snp.filter.vcf
gatk --java-options "-Xmx70G" SelectVariants -R $ref -V ref.snp.filter.vcf --exclude-filtered -O ref.snp.final.vcf.gz

# 04. phasing
for i in {1..19}; do
cat << EOF > chr_${i}.shapeit.sh 
#!/bin/bash
#SBATCH --job-name=shapeit_chr${i}
#SBATCH --time=100:00:00 # hh:mm:ss
#SBATCH --partition=cpu
#SBATCH -N 1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=20G
#SBATCH --output=%j.out
#SBATCH --error=%j.err

# 011. split chrs
vcftools --gzvcf ref.snp.final.vcf.gz --chr ${i}  --maf 0.05 --recode --recode-INFO-all --out chr${i}
# 012. beagle
## change .(the missing data) to ./.(the missing genotype)
cat chr${i}.recode.vcf | perl -pe "s/\s\.:/\t.\/.:/g" > chr${i}.recode2.vcf
java -Xmx150G -jar /share/home/zhou3lab/suqianwei/tools/beagle.22Jul22.46e.jar gt=chr${i}.recode2.vcf out=chr${i}_snp.IM nthreads=4

gunzip chr${i}_snp.IM.vcf.gz
# only select snp with one nucletide in alt position
awk '\$1~/#/ || \$5!~/,/' chr${i}_snp.IM.vcf > chr${i}_snp.IM.bi.vcf
# 013. shapeit phase
/share/home/zhou3lab/suqianwei/tools/shapeit.v2.904.3.10.0-693.11.6.el7.x86_64/bin/shapeit -V chr${i}_snp.IM.bi.vcf --window 0.5 -O chr${i}_snp.IM.Ph --thread 4
EOF
sbatch chr_${i}.shapeit.sh
done

