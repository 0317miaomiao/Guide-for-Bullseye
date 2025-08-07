#!/bin/bash
#SBATCH --job-name=dedup_bam_parallel
#SBATCH -p kshcnormal
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=32
#SBATCH --time=7-00:00:00
#SBATCH --output=dedup_%j.out
#SBATCH --error=dedup_%j.err

source activate bullseye_env

INPUT_DIR="bam_output"
OUTPUT_DIR="bam_sort"

mkdir -p "$OUTPUT_DIR"

# 定义处理每个SRR目录的函数
process_sample() {
    SAMPLE_ID="$1"
    SAMPLE_DIR="$INPUT_DIR/$SAMPLE_ID"
    RAW_BAM="$SAMPLE_DIR/${SAMPLE_ID}_Aligned.sortedByCoord.out.bam"

    if [[ ! -f "$RAW_BAM" ]]; then
        echo "❌ BAM file not found: $RAW_BAM"
        return
    fi

    echo "🔄 Processing $SAMPLE_ID..."

    OUT_DIR="$OUTPUT_DIR/$SAMPLE_ID"
    mkdir -p "$OUT_DIR"

    # Step 1: name sort
    samtools sort -n -o "$OUT_DIR/${SAMPLE_ID}_namesorted.bam" "$RAW_BAM"

    # Step 2: fixmate
    samtools fixmate -m "$OUT_DIR/${SAMPLE_ID}_namesorted.bam" "$OUT_DIR/${SAMPLE_ID}_fixmate.bam"

    # Step 3: coordinate sort
    samtools sort -o "$OUT_DIR/${SAMPLE_ID}_coord_sorted.bam" "$OUT_DIR/${SAMPLE_ID}_fixmate.bam"

    # Step 4: mark duplicates
    samtools markdup -r "$OUT_DIR/${SAMPLE_ID}_coord_sorted.bam" "$OUT_DIR/${SAMPLE_ID}_dedup.bam"

    # Step 5: index the final deduplicated BAM
    samtools index "$OUT_DIR/${SAMPLE_ID}_dedup.bam"

    # Step 6: cleanup intermediate files
    rm "$OUT_DIR/${SAMPLE_ID}_namesorted.bam"
    rm "$OUT_DIR/${SAMPLE_ID}_fixmate.bam"
    rm "$OUT_DIR/${SAMPLE_ID}_coord_sorted.bam"

    echo "✅ Done with $SAMPLE_ID"
}

export -f process_sample
export INPUT_DIR OUTPUT_DIR

# 获取所有 SRR 子目录名（如 SRR15268477）写入临时文件
find "$INPUT_DIR" -maxdepth 1 -type d -name 'SRR*' | xargs -n 1 basename > srr_list_dedup.txt

# 并行执行，每次最多 4 个任务
cat srr_list_dedup.txt | xargs -P 8 -n 1 -I {} bash -c 'process_sample "$@"' _ {}