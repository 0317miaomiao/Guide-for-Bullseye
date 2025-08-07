#!/bin/bash
#SBATCH --job-name=bullseye_merge_2
#SBATCH -p kshcnormal
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=32
#SBATCH --time=7-00:00:00
#SBATCH --output=bullseye_merge_2_%j.out
#SBATCH --error=bullseye_merge_2_%j.err

# 激活环境
source activate bullseye_env

perl summarize_sites.pl --minRep 324 result_merge/*.bed > result_merge/result.bed