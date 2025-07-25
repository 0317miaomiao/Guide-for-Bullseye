# Guide-for-Bullseye

This guide mainly describes how to use Bullseye to perform the overall workflow of single-cell methylation quantification. It first briefly introduces the main installation steps and core functions of Bullseye, and then demonstrates the complete process using the HEK293T dataset as an example. The full code and corresponding comments can be found in the repository.

## How to Install and Use Bullseye
### 1. GitHub Repository
```bash
git clone https://github.com/mflamand/Bullseye.git
```

### 2. Environment Setup
Create a new conda environment:

```bash
conda create -n bullseye_env
```
Install required packages:

```bash
# Install Perl
conda install -c conda-forge perl=5.26

# Install samtools
conda install -c bioconda samtools=1.12

# Install htslib
conda install -c bioconda htslib

# Install bedtools
conda install -c bioconda bedtools

# Install Bio::DB::Fasta
conda install -c bioconda perl-bioperl

# Install required Perl modules
cpanm MCE
cpanm Array::IntSpan
cpanm Text::NSP
```
**Note: It is strongly recommended to use the specified software versions to avoid compatibility issues.**

### 3. Installation Check

Check software versions:
```bash
perl -v            # Should show Perl version > 5.26
samtools --version # Should display samtools version info
bedtools --version # Should display bedtools version info
tabix --version    # Should display tabix version info
```

Check Perl modules:

```bash
perl -MMCE -e1
perl -MMath::CDF -e1
perl -MBio::DB::Fasta -e1
perl -MText::NSP -e1
perl -MArray::IntSpan -e1
perl -MXML::Parser -e1
```

If the module is installed correctly, nothing will be printed (successful silent exit).
If not installed, an error like ```Can't locate XYZ.pm in @INC ...``` will be shown.


### 4. Usage Guide
#### 4.1 Quantification
```bash
perl parseBAM.pl --input file.bam --output output.matrix --cpu 4 --minCoverage 10 --removeDuplicates
```

The input ```.bam``` file must be sorted and accompanied by its ```.bai ``` index.

The output will include output.matrix.gz and its index file output.matrix.gz.tbi.

```bash
tabix -s 1 -b 2 -e 2 output.matrix.gz
```

#### 4.2 RNA Editing Site Detection
```bash
perl Find_edit_site.pl --annotationFile $annotation_file \
         --EditedMatrix dart.matrix.gz \
         --controlMatrix control.matrix.gz \
         --minEdit 5 \              # minimal editing rate
         --maxEdit 90 \             # maximal editing rate
         --editFoldThreshold 1.5 \  # minimum editing fold change over control
         --MinEditSites 3 \         # minimum number of mutations required for detection
         --cpu 4 \
         --outfile output.bed \
         --fallback genome.fasta \  # reference genome, used when a site lacks coverage in control
         --verbose
```
```annotation_file``` should be a refFlat file.

Please refer to the [GitHub repository](https://github.com/mflamand/Bullseye) for instructions on creating this file.

Currently supported species: Human and Mouse.

## HEK293T Download and Processing Workflow

### 1. Download Related Data

Mainly use single-cell data, including iDART_SMART-seq APOBEC1-YTH Cell and iDART_SMART-seq APOBEC1-YTHmut Cell. Separate them into two folders based on whether they are mutants (mut), then use download.sh and download_mut.sh scripts to download separately. Note to modify the parallel count according to your server configuration.

### 2. Data Processing

#### 2.1 Use Flexbar for Adapter Trimming and Quality Cutting (flexbar.sh and flexbar_mut.sh)

```bash
flexbar \
  -r /path/to/..._1.fastq.gz \
  -p /path/to/..._2.fastq.gz \
  -t trimmed/sample_trimmed \
  --adapter-preset Nextera \
  -ap ON \
  --adapter-trim-end RIGHT \
  --adapter-min-overlap 7 \
  -ae 0.1 \
  --qtrim-format sanger \
  -q TAIL \
  -qtrim-threshold 20 \
  -m 36 \
  -n 8 \
  --zip-output GZ
```

**Parameter Explanations:**
- `-r /path/to/..._1.fastq.gz`: Input file 1 (Read 1), path to the first end FASTQ file of paired-end sequencing
- `-p /path/to/..._2.fastq.gz`: Input file 2 (Read 2), path to the second end FASTQ file of paired-end sequencing
- `-t trimmed/sample_trimmed`: Output file prefix, processed files will be named with this prefix, e.g., trimmed/sample_trimmed_1.fastq.gz
- `--adapter-preset Nextera`: Use built-in Nextera adapter preset, automatically loads Illumina Nextera library construction adapter sequences
- `-ap ON`: Enable paired-end adapter overlap detection, use paired information for more accurate adapter removal
- `--adapter-trim-end RIGHT`: Only remove adapters from the right end (3' end) of reads, as Illumina adapters are usually connected at sequence ends
- `--adapter-min-overlap 7`: Perform trimming only when adapter and read match with at least 7 base overlap, avoiding false trimming
- `-ae 0.1`: Allow 10% mismatches (error rate 0.1), permit certain mismatches when adapter matches read
- `--qtrim-format sanger`: Specify FASTQ quality value format as Sanger (Phred+33), common format for modern Illumina sequencing
- `-q TAIL`: Perform quality trimming from read tail (low-quality bases will be trimmed)
- `-qtrim-threshold 20`: Quality threshold of 20, bases below this Phred score will be trimmed from the tail
- `-m 36`: Minimum retained length 36 bp, discard reads shorter than 36 bases after trimming
- `-n 8`: Use 8 threads for parallel processing to speed up
- `--zip-output GZ`: Compress output files to gzip format to reduce file size

### 3. Use STAR to Generate BAM Files

Generate genome index (star_1.sh):
```bash
STAR --runThreadN 8 \
     --runMode genomeGenerate \
     --genomeDir genome_index \
     --genomeFastaFiles /path/to/...human.fa \
     --sjdbGTFfile /path/to/...Homo_sapiens.GRCh38.112.gtf \
     --sjdbOverhang 150
```

This generates genome_index for the next alignment step. The `--sjdbOverhang` parameter is determined by the read length.  
Specifically, it should be set to **read length minus 1**。

STAR alignment (star_2.sh and star_2mut.sh):
```bash
STAR --runThreadN 8 \
     --genomeDir genome_index \
     --readFilesIn "$R1" "$R2" \
     --readFilesCommand zcat \
     --outFileNamePrefix "$OUTPUT_SUBDIR/${SRR_ID}_" \
     --outSAMtype BAM SortedByCoordinate
```
### 4. BAM File Processing and Deduplication using SAMtools

After generating BAM files with STAR, perform deduplication using SAMtools (1.11) fixmate and markdup:

```bash
# Step 1: Name sort the BAM file
samtools sort -n -o path/to/bam_output/${SAMPLE_ID}_namesorted.bam path/to/star_output/${SAMPLE_ID}_Aligned.sortedByCoord.out.bam

# Step 2: Fix mate information
samtools fixmate -m path/to/bam_output/${SAMPLE_ID}_namesorted.bam path/to/bam_output/${SAMPLE_ID}_fixmate.bam

# Step 3: Coordinate sort the fixed BAM file
samtools sort -o path/to/bam_output/${SAMPLE_ID}_coord_sorted.bam path/to/bam_output/${SAMPLE_ID}_fixmate.bam

# Step 4: Mark and remove duplicates
samtools markdup -r path/to/bam_output/${SAMPLE_ID}_coord_sorted.bam path/to/bam_output/${SAMPLE_ID}_dedup.bam

# Step 5: Index the final deduplicated BAM file
samtools index path/to/bam_output/${SAMPLE_ID}_dedup.bam
```

**Process Explanation:**
1. **Name sort**: Sort BAM file by read names to prepare for fixmate
2. **Fixmate**: Fix mate pair information and add mate score tags
3. **Coordinate sort**: Sort BAM file by genomic coordinates for duplicate marking
4. **Mark duplicates**: Identify and remove PCR/optical duplicates using `-r` flag
5. **Index**: Create BAM index file for efficient random access

The final output `${SAMPLE_ID}_dedup.bam` is the processed BAM file ready for downstream analysis.













