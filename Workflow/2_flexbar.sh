#!/bin/bash
#SBATCH --job-name=flexbar_parallel
#SBATCH -p kshcnormal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16           # 分配16核
#SBATCH --time=7-00:00:00
#SBATCH --output=batch_flexbar%j.out
#SBATCH --error=batch_flexbar%j.err

source activate bullseye_pre

FASTQ_DIR="APOBEC1-YTH_fastq"
OUTPUT_DIR="trimmed"

mkdir -p "$OUTPUT_DIR"

download_srr() {
    srr_id="$1"
    srr_path="${FASTQ_DIR}/${srr_id}"
    fq1="${srr_path}/${srr_id}_1.fastq.gz"
    fq2="${srr_path}/${srr_id}_2.fastq.gz"

    if [[ ! -f "$fq1" || ! -f "$fq2" ]]; then
        echo "[$srr_id] Missing FASTQ files, skipping."
        return
    fi

    outdir="${OUTPUT_DIR}/${srr_id}"
    mkdir -p "$outdir"
    outprefix="${outdir}/${srr_id}_trimmed"

    echo "[$srr_id] Processing..."

    flexbar \
      -r "$fq1" \
      -p "$fq2" \
      -t "$outprefix" \
      --adapter-preset Nextera \
      -ap ON \
      --adapter-trim-end RIGHT \
      --adapter-min-overlap 7 \
      -ae 0.1 \
      --qtrim-format sanger \
      -q TAIL \
      -qtrim-threshold 20 \
      -m 36 \
      -n 2 \
      --zip-output GZ

    if [[ $? -ne 0 ]]; then
        echo "[$srr_id] Flexbar failed. Removing output directory."
        rm -rf "$outdir"
    else
        echo "[$srr_id] Finished successfully."
    fi
}

export -f download_srr
export FASTQ_DIR OUTPUT_DIR

# 找所有SRR目录名称，写入临时文件
find "$FASTQ_DIR" -maxdepth 1 -type d -name 'SRR*' | xargs -n 1 basename > srr_list.txt

cat srr_list.txt | xargs -P 4 -n 1 -I {} bash -c 'download_srr "$@"' _ {}
