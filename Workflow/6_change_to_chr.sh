#!/bin/bash
#SBATCH --job-name=change_to_chr
#SBATCH -p kshcnormal
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=32
#SBATCH --time=72:00:00
#SBATCH -o change_to_chr_%j.out
#SBATCH -e change_to_chr_%j.err

source activate bullseye_env

# 定义输入和输出文件夹
input_dir="bull_mutmatrix"
output_dir="bull_mutmatrix_chr"

# 创建输出文件夹
mkdir -p "$output_dir"

# 遍历输入文件夹中的所有子文件夹
for subdir in "$input_dir"/*/; do
    # 获取子文件夹名称（如SRR15268425）
    sample_name=$(basename "$subdir")
    
    # 定义输入和输出文件路径
    input_gz="${subdir}${sample_name}.matrix.gz"
    output_gz="${output_dir}/${sample_name}.matrix.gz"
    
    # 检查输入文件是否存在
    if [ -f "$input_gz" ]; then
        echo "正在处理 $input_gz"
        
        # 解压.gz文件，添加chr前缀，然后重新压缩
        gunzip -c "$input_gz" | awk '{ $1 = "chr"$1; print }' OFS="\t" | bgzip -c > "$output_gz"
        
        # 为新的.gz文件创建.tbi索引
        tabix -s 1 -b 2 -e 2 "$output_gz"
        
        echo "处理完成: ${sample_name}.matrix.gz 和 ${sample_name}.matrix.gz.tbi 已生成"
    else
        echo "警告: $input_gz 不存在，跳过"
    fi
done

echo "所有文件处理完成，结果保存在 $output_dir 文件夹中"