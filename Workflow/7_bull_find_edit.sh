#!/bin/bash
#SBATCH --job-name=bullseye_parallel_grouped
#SBATCH -p kshcnormal
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=32
#SBATCH --time=7-00:00:00
#SBATCH --output=batch_edit_%j.out
#SBATCH --error=batch_edit_%j.err

# 激活环境
source activate bullseye_env

# 设置变量
WT_DIR="bull_matrix_chr"
MUT_DIR="bull_mutmatrix_chr"
OUT_DIR="result"
SCRIPT="Find_edit_site.pl"
ANNOTATION="Prepare/human.refFlat"

mkdir -p "$OUT_DIR"

# 定义函数：一个 WT 文件对应多个 MUT 文件
process_group() {
    WT_FILE="$1"
    WT_NAME=$(basename "$WT_FILE" .matrix.gz)
    SUB_DIR="$OUT_DIR/$WT_NAME"
    mkdir -p "$SUB_DIR"

    count=1
    for MUT_FILE in "$MUT_DIR"/*.matrix.gz; do
        OUTFILE="$SUB_DIR/result${count}.bed"

        perl "$SCRIPT" \
            --annotationFile "$ANNOTATION" \
            --EditedMatrix "$WT_FILE" \
            --controlMatrix "$MUT_FILE" \
            --minEdit 5 \
            --maxEdit 90 \
            --editFoldThreshold 1.5 \
            --MinEditSites 2 \
            --EditedMinCoverage 5 \
            --ControlMinCoverage 5 \
            --cpu 8 \
            --outfile "$OUTFILE" \
            --verbose

        count=$((count + 1))
    done
}

export -f process_group
export SCRIPT ANNOTATION MUT_DIR OUT_DIR

# 使用 xargs 并行处理每个 WT 对应的全部 MUT 文件
find "$WT_DIR" -name "*.matrix.gz" | xargs -P 6 -n 1 -I {} bash -c 'process_group "$@"' _ {}