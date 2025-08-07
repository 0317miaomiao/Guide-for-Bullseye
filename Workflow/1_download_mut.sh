#!/bin/bash
#SBATCH --job-name=download_sra
#SBATCH -p kshcnormal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=7-00:00:00
#SBATCH --output=mut_%j.out
#SBATCH --error=mut_%j.err


# 设置输入文件和输出目录
SRR_LIST="APOBEC1-YTHmut_SRR_IDs.txt"
OUT_DIR="APOBEC1-YTHmut_fastq"

# 创建主输出目录
mkdir -p "$OUT_DIR"

# 定义每个 SRR 的下载逻辑（作为函数，供 xargs 调用）
download_srr() {
    SRR_ID=$1
    echo "🔄 Processing $SRR_ID..."

    # 构建子目录
    SUB_DIR="${OUT_DIR}/${SRR_ID}"
    mkdir -p "$SUB_DIR"

    # 构造 FTP 路径（取最后两位并补0）
    PREFIX=${SRR_ID:0:6}
    LAST2=${SRR_ID: -2}
    SUBDIR="0${LAST2}"
    FTP_PATH="https://ftp.sra.ebi.ac.uk/vol1/fastq/${PREFIX}/${SUBDIR}/${SRR_ID}"

    # 下载 _1.fastq.gz
    FILE1="${SRR_ID}_1.fastq.gz"
    URL1="${FTP_PATH}/${FILE1}"
    echo "  ⬇️ Downloading $FILE1..."
    wget -q -c "$URL1" -O "${SUB_DIR}/${FILE1}"

    if [ ! -s "${SUB_DIR}/${FILE1}" ]; then
        echo "  ❌ Failed to download ${FILE1}. Skipping."
        return
    fi

    # 检查 _2.fastq.gz 是否存在
    FILE2="${SRR_ID}_2.fastq.gz"
    URL2="${FTP_PATH}/${FILE2}"
    echo "  🔍 Checking for $FILE2..."

    if wget --spider -q "$URL2"; then
        echo "  ⬇️ Downloading $FILE2..."
        wget -q -c "$URL2" -O "${SUB_DIR}/${FILE2}"
    else
        echo "  ⚠️  $FILE2 not found. Single-end only."
    fi

    echo "✅ Finished $SRR_ID."
    echo "----------------------------"
}

export -f download_srr
export OUT_DIR

# 控制并行数（可根据服务器核数调整，比如你有 64 核，建议用 16 或 32）
cat "$SRR_LIST" | xargs -P 16 -n 1 -I {} bash -c 'download_srr "$@"' _ {}

