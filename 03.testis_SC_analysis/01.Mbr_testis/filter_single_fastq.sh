#!/bin/bash
#SBATCH -J f2
#SBATCH -p cpu
#SBATCH -N 1
#SBATCH --mem=10G
#SBATCH --ntasks-per-node=4
#SBATCH --output=%J.out
#SBATCH --error=%J.err


# 用法: bash filter_single_fastq.sh <input_fastq.gz> <output_fastq.gz> <remove_id_list.txt>

set -e

INPUT_FASTQ="$1"
OUTPUT_FASTQ="$2"
REMOVE_LIST="$3"

# 自动识别是 R1 还是 R2
if [[ "$INPUT_FASTQ" =~ \/1\.fq\.gz$ ]] || [[ "$INPUT_FASTQ" =~ _1\.fq\.gz$ ]]; then
    SUFFIX="/1"
elif [[ "$INPUT_FASTQ" =~ \/2\.fq\.gz$ ]] || [[ "$INPUT_FASTQ" =~ _2\.fq\.gz$ ]]; then
    SUFFIX="/2"
else
    echo "无法判断是 R1 还是 R2，请检查文件名: $INPUT_FASTQ"
    exit 1
fi

echo "[INFO] Processing $INPUT_FASTQ as read $SUFFIX"

# 统计原始 reads 数
TOTAL=$(zcat "$INPUT_FASTQ" | awk 'NR % 4 == 1' | wc -l)

# 过滤 reads
seqkit grep -v -f "$REMOVE_LIST" -j 4 "$INPUT_FASTQ" -o "$OUTPUT_FASTQ"

# 统计过滤后 reads 数
REMAINING=$(zcat "$OUTPUT_FASTQ" | awk 'NR % 4 == 1' | wc -l)
REMOVED=$((TOTAL - REMAINING))

echo "[STAT] Total reads:     $TOTAL"
echo "[STAT] Removed reads:   $REMOVED"
echo "[STAT] Remaining reads: $REMAINING"

# 清理
rm "$TMP_REMOVE_LIST"

