#!/bin/bash
#SBATCH --job-name=run_parseBAMmut
#SBATCH -p kshcnormal
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=32
#SBATCH --time=72:00:00
#SBATCH -o parseBAMmut_%j.out
#SBATCH -e parseBAMmut_%j.err

# 配置环境变量
BAM_DIR="bam_mutsort"
OUTPUT_DIR="bull_mutmatrix"
SCRIPT_PATH="parseBAM.pl"

# 创建输出目录
mkdir -p "$OUTPUT_DIR"

# 激活 conda 环境
source activate bullseye_env

# 定义函数（对单个SRR运行）
process_bam() {
    SRR_ID=$1
    INPUT_BAM="${BAM_DIR}/${SRR_ID}/${SRR_ID}_dedup.bam"
    OUTPUT_FOLDER="${OUTPUT_DIR}/${SRR_ID}"
    OUTPUT_MATRIX="${OUTPUT_FOLDER}/${SRR_ID}.matrix"

    echo "🔧 Processing $SRR_ID"

    # 如果 BAM 文件存在才处理
    if [[ -f "$INPUT_BAM" ]]; then
        mkdir -p "$OUTPUT_FOLDER"
        perl "$SCRIPT_PATH" \
            --input "$INPUT_BAM" \
            --output "$OUTPUT_MATRIX" \
            --cpu 8 \
            --minCoverage 10 \
            --removeDuplicates
    else
        echo "❌ BAM file not found for $SRR_ID: $INPUT_BAM"
    fi
}

export -f process_bam
export BAM_DIR OUTPUT_DIR SCRIPT_PATH

# 获取所有 SRR 子目录名
find "$BAM_DIR" -maxdepth 1 -type d -name 'SRR*' | xargs -n 1 basename > bull_mutmatrix.txt

cat bull_mutmatrix.txt | xargs -P 6 -n 1 -I {} bash -c 'process_bam "$@"' _ {}

rm bull_mutmatrix.txt
