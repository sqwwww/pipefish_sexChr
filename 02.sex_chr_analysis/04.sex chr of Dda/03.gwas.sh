#!/bin/bash
#SBATCH --job-name=gwas
#SBATCH --time=70:00:00 # hh:mm:ss
#SBATCH --partition=cpu
#SBATCH -N 1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=50G
#SBATCH --output=%j.out
#SBATCH --error=%j.err

for i in {1..21}; do cat chr${i}_snp.IM.Ph.haps >>snp.haps; done
# change snp.sample to correct format for getsnpplink
cp chr1_snp.IM.Ph.sample snp.sample
#awk '{if(NR<3)print $0; if (NR>2){split($1,a,"-");print a[3]"_"$1,a[3]"_"$1,0}}' snp.sample | sed 's/^\s*$//' > snp.sample.correct

# convert format
perl /share/home/zhou3lab/suqianwei/tools/getsnpplink.pl snp.sample snp.haps snp.map snp.ped
plink --file snp --out snp --make-bed --allow-extra-chr --allow-no-sex

# generate tped/tfam
plink --bfile snp --autosome-num 21 --allow-no-sex --recode 12 transpose --output-missing-genotype 0 --out snp
# generate pca data
plink --bfile snp --autosome-num 21 --allow-no-sex --pca --out snp

# generate phenotype file
awk '{print $1,$2,$6}' snp.tfam > snp.pheno
# replace character in chr to label in tped and tfam
#sed 's/_hap1//g' snp.tped|sed 's/chr//g' > snp.correct.tped
#cp snp.tfam snp.correct.tfam
# run emmax and use snp.hBN as final result.
/share/home/zhou3lab/suqianwei/tools/emmax-kin-intel64 -v -s -d 10 snp
/share/home/zhou3lab/suqianwei/tools/emmax-kin-intel64 -v -d 10 snp
/share/home/zhou3lab/suqianwei/tools/emmax-intel64 -v -d 10 -t snp -p snp.pheno -k snp.aIBS.kinf -o snp.aIBS
/share/home/zhou3lab/suqianwei/tools/emmax-intel64 -v -d 10 -t snp -p snp.pheno -k snp.aBN.kinf -o snp.aBN

awk '{split($1,a,":");print $1,a[1],a[2],$4}' snp.aBN.ps > dataForManhattan.aBN.txt
awk '{split($1,a,":");print $1,a[1],a[2],$4}' snp.aIBS.ps > dataForManhattan.aIBS.txt
sbatch plotMH.sh
