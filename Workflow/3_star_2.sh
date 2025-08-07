#!/bin/bash
#SBATCH --job-name=star_batch
#SBATCH -p kshcnormal
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=32
#SBATCH --time=7-00:00:00
#SBATCH --output=star2_%j.out
#SBATCH --error=star2_%j.err

source activate bullseye_pre

# 设置输入输出目录
FASTQ_DIR="trimmed"
OUTPUT_DIR="bam_output"

# 定义对每个 SRR 样本的处理函数
download_srr() {
    SRR_ID="$1"
    INPUT_DIR="$FASTQ_DIR/$SRR_ID"
    OUTPUT_SUBDIR="$OUTPUT_DIR/$SRR_ID"

    # 如果输出文件夹已经存在，说明已经跑过，跳过
    if [ -d "$OUTPUT_SUBDIR" ]; then
        echo "✅ 已存在输出目录 $OUTPUT_SUBDIR，跳过 $SRR_ID"
        return
    fi

    echo "🔄 开始处理 $SRR_ID..."

    # 创建输出目录
    mkdir -p "$OUTPUT_SUBDIR"

    # 获取 FASTQ 文件路径
    R1="$INPUT_DIR/${SRR_ID}_trimmed_1.fastq.gz"
    R2="$INPUT_DIR/${SRR_ID}_trimmed_2.fastq.gz"

    # 检查输入文件是否存在
    if [[ ! -f "$R1" || ! -f "$R2" ]]; then
        echo "⚠️ 缺少 $R1 或 $R2，跳过 $SRR_ID"
        return
    fi

    # STAR 比对（根据需要替换参数）
    STAR --runThreadN 8 \
         --genomeDir genome_index \
         --readFilesIn "$R1" "$R2" \
         --readFilesCommand zcat \
         --outFileNamePrefix "$OUTPUT_SUBDIR/${SRR_ID}_" \
         --outSAMtype BAM SortedByCoordinate

    echo "✅ 完成 $SRR_ID"
}

export -f download_srr
export FASTQ_DIR OUTPUT_DIR

# 获取所有 SRR 子目录名（如 SRR15268477）写入临时文件
find "$FASTQ_DIR" -maxdepth 1 -type d -name 'SRR*' | xargs -n 1 basename > srr_list_star.txt

# 并行执行，每次最多 4 个任务
cat srr_list_star.txt | xargs -P 3 -n 1 -I {} bash -c 'download_srr "$@"' _ {}