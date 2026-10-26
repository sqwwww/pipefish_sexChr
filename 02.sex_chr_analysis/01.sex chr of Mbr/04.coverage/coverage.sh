#!/bin/bash
#SBATCH --job-name=coverage
#SBATCH --partition=cpu
#SBATCH -N 1
#SBATCH --mem=20G
#SBATCH --ntasks-per-node=4
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH --time=50:00:00

# Calculate male and female sequencing depth in 10-kb genomic windows.
# Reads from 10 males and 10 females are pooled separately.
# Repeat regions are excluded before calculating window-level statistics.
#
# Requirements: SAMtools and BEDTools.
# Input BAM files must be coordinate-sorted and use the same reference.
# The repeat annotation must be in BED format and use chromosome names
# consistent with the BAM files and reference index.
#
# These sections can also be run separately if intermediate files exist.

# Stop if a command fails.
set -e
set -o pipefail

# Input files: adjust paths as needed.
genome_fai=../HJ.genome.masked.chr.fa.fai
repeat_bed=HJ.genome.chr.fasta.out.gff.bed
loc_bed=HJ.genome.masked.chr.fa.10k.bed


# 1. Generate non-overlapping 10-kb windows.
# The final window of each chromosome may be shorter than 10 kb.
bedtools makewindows \
    -g "$genome_fai" \
    -w 10000 > "$loc_bed"


# 2. Merge individual BAM files by sex and index the pooled BAM files.
samtools merge male.bam \
    ../M1.bam ../M2.bam ../M3.bam ../M4.bam ../M5.bam \
    ../M6.bam ../M7.bam ../M8.bam ../M9.bam ../M10.bam

samtools index male.bam

samtools merge female.bam \
    ../F1.bam ../F2.bam ../F3.bam ../F4.bam ../F5.bam \
    ../F6.bam ../F7.bam ../F8.bam ../F9.bam ../F10.bam

samtools index female.bam


# 3. Calculate per-base sequencing depth, including zero-depth positions.
samtools depth -aa male.bam > male_samtools.depth
samtools depth -aa female.bam > female_samtools.depth


# 4. Convert 1-based positions to 0-based, half-open BED intervals.
# Sort both the windows and depth files using the same chromosome order.
export LC_ALL=C

sort -k1,1 -k2,2n "$loc_bed" > "${loc_bed}.tmp"
mv "${loc_bed}.tmp" "$loc_bed"

awk 'BEGIN {OFS="\t"} {print $1, $2-1, $2, $3}' \
    male_samtools.depth |
    sort -k1,1 -k2,2n > male_samtools.depth.bed

awk 'BEGIN {OFS="\t"} {print $1, $2-1, $2, $3}' \
    female_samtools.depth |
    sort -k1,1 -k2,2n > female_samtools.depth.bed


# 5. Remove positions overlapping annotated repeat regions.
bedtools subtract \
    -a male_samtools.depth.bed \
    -b "$repeat_bed" \
    > male_samtools_subtractRepeat.depth.bed

bedtools subtract \
    -a female_samtools.depth.bed \
    -b "$repeat_bed" \
    > female_samtools_subtractRepeat.depth.bed


# 6. Summarize depth across non-repeat positions in each 10-kb window.
# Output columns:
# 1-3: chromosome, window start, window end
# 4: sum of per-base depths
# 5: number of non-repeat positions with depth records (including zeros)
# 6: median depth
# 7: mean depth
# Windows without eligible positions are reported as ".".
bedtools map \
    -a "$loc_bed" \
    -b male_samtools_subtractRepeat.depth.bed \
    -c 4 -o sum,count,median,mean \
    > male_samtools.10kwindow_Srepeat.depth

bedtools map \
    -a "$loc_bed" \
    -b female_samtools_subtractRepeat.depth.bed \
    -c 4 -o sum,count,median,mean \
    > female_samtools.10kwindow_Srepeat.depth


# 7. Optional: extract Chr1 alignments for downstream analyses.
# Uncomment to run. Use the chromosome name present in the BAM header.
# Here, chromosome 1 is named "1".
#
# samtools view -@ 4 -b -o female.chr1.bam female.bam 1
# samtools index -@ 4 female.chr1.bam
#
# samtools view -@ 4 -b -o male.chr1.bam male.bam 1
# samtools index -@ 4 male.chr1.bam


# 8. Optional: count non-repeat positions with pooled depth >= 3.
# Uncomment to run. The fourth output column is the number of positions.
# A value of 0 is reported for windows with no qualifying positions.
#
# awk '$4 >= 3' male_samtools_subtractRepeat.depth.bed \
#     > male_samtools.SR.cov3.depth.bed
#
# bedtools map \
#     -a "$loc_bed" \
#     -b male_samtools.SR.cov3.depth.bed \
#     -c 4 -o count -null 0 \
#     > male_samtools.10kwindow.SR.cov3.coveredSites.bed
#
# awk '$4 >= 3' female_samtools_subtractRepeat.depth.bed \
#     > female_samtools.SR.cov3.depth.bed
#
# bedtools map \
#     -a "$loc_bed" \
#     -b female_samtools.SR.cov3.depth.bed \
#     -c 4 -o count -null 0 \
#     > female_samtools.10kwindow.SR.cov3.coveredSites.bed