#!/bin/bash
#SBATCH --job-name=bullseye_merge_1
#SBATCH -p kshcnormal
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=32
#SBATCH --time=7-00:00:00
#SBATCH --output=batch_merge_1_%j.out
#SBATCH --error=batch_merge_1_%j.err

# 激活环境
source activate bullseye_env

# 路径变量
INPUT_DIR="result_new"
OUTPUT_DIR="result_merge"
SUMMARIZE_SCRIPT="summarize_sites.pl"

# 创建输出目录
mkdir -p "$OUTPUT_DIR"

# 定义函数用于处理一个子文件夹
process_folder() {
    FOLDER_NAME="$1"
    BED_DIR="$INPUT_DIR/$FOLDER_NAME"
    OUTPUT_FILE="$OUTPUT_DIR/${FOLDER_NAME}.bed"

    # 构建所有 .bed 文件路径
    BED_FILES=("$BED_DIR"/*.bed)

    # 如果没有 .bed 文件就跳过
    if [ ${#BED_FILES[@]} -eq 0 ]; then
        echo "No BED files found in $BED_DIR"
        return
    fi

    # 运行合并命令
    perl "$SUMMARIZE_SCRIPT" --repOnly "${BED_FILES[@]}" > "$OUTPUT_FILE"
}

# 导出函数和变量供 xargs 使用
export -f process_folder
export INPUT_DIR OUTPUT_DIR SUMMARIZE_SCRIPT

# 获取所有 SRR 文件夹名并并行处理
find "$INPUT_DIR" -mindepth 1 -maxdepth 1 -type d -printf "%f\n" | \
    xargs -n 1 -P 12 -I {} bash -c 'process_folder "$@"' _ {}