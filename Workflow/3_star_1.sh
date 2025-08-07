#!/bin/bash
#SBATCH --job-name=star_batch
#SBATCH -p kshcnormal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=7-00:00:00
#SBATCH --output=star1_%j.out
#SBATCH --error=star1_%j.err

source activate bullseye_pre

STAR --runThreadN 8 \
     --runMode genomeGenerate \
     --genomeDir genome_index \
     --genomeFastaFiles Prepare/index/merged_chromosomes.fa \
     --sjdbGTFfile Prepare/human_gtf/Homo_sapiens.GRCh38.112.gtf \
     --sjdbOverhang 150
