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

### 5. Generate Count Matrix Using Bullseye

Use the following command to generate a count matrix through Bullseye:

```bash
perl parseBAM.pl \
    --input "$INPUT_BAM" \
    --output "$OUTPUT_MATRIX" \
    --cpu 8 \
    --minCoverage 10 \
    --removeDuplicates
```

#### Parameter Description

- `--input "$INPUT_BAM"`: Specify the input BAM file path. This BAM file should be aligned and sorted sequencing data.

- `--output "$OUTPUT_MATRIX"`: Specify the output count matrix file path. The generated matrix will contain expression counts for each gene in each cell.

- `--cpu 8`: Set the number of CPU cores for parallel processing. Set to 8 cores here, can be adjusted according to server configuration to improve processing speed.

- `--minCoverage 10`: Set the minimum coverage threshold. Only sites with coverage of 10 or above will be included in the analysis to filter low-quality data.

- `--removeDuplicates`: Flag to remove PCR duplicate sequences. Enabling this option removes duplicate reads generated by PCR amplification, improving data quality and analysis accuracy.

### 6. Add chr Prefix

Since the standard Bullseye code requires the first column of the matrix generated in step 5 to be in "chr1" format (i.e., with "chr" prefix), and our annotation file may not contain this prefix, modifications are needed, otherwise Bullseye cannot be used correctly.

#### Processing Steps

Use the following commands to add chr prefix to the matrix file:

```bash
# Decompress .gz file, add chr prefix, then recompress
gunzip -c "$input_gz" | awk '{ $1 = "chr"$1; print }' OFS="\t" | bgzip -c > "$output_gz"

# Create .tbi index for the new .gz file
tabix -s 1 -b 2 -e 2 "$output_gz"
```

**Important**: After adding the chr prefix, you need to regenerate the index file (.tbi) to ensure the file can be properly indexed and accessed.

#### Command Description

- `gunzip -c "$input_gz"`: Decompress the input .gz file and output to standard output
- `awk '{ $1 = "chr"$1; print }' OFS="\t"`: Add "chr" prefix to the first column and use tab as output separator
- `bgzip -c > "$output_gz"`: Recompress using bgzip and save to output file
- `tabix -s 1 -b 2 -e 2 "$output_gz"`: Create index for compressed file, specifying column 1 as sequence name, column 2 as start position, and column 2 as end position

### 7. Generate Final Statistical Results Using Bullseye

Compare each wild-type file with all mutant files to calculate final quantitative results and identify editing sites:

```bash
perl Find_edit_site.pl \
    --annotationFile human.refFlat \
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
```

#### Parameter Description

- `--annotationFile human.refFlat`: Specify the gene annotation file to determine the position and function of editing sites within genes.

- `--EditedMatrix "$WT_FILE"`: Specify the count matrix file path for the edited group (wild-type).

- `--controlMatrix "$MUT_FILE"`: Specify the count matrix file path for the control group (mutant).

- `--minEdit 5`: Set minimum editing rate threshold to 5%. Sites below this threshold will not be considered as editing sites.

- `--maxEdit 90`: Set maximum editing rate threshold to 90%. Sites above this threshold may have genotype variations rather than RNA editing.

- `--editFoldThreshold 1.5`: Set editing fold change threshold. The ratio of editing rates between edited and control groups must be greater than 1.5-fold to be considered a significant editing site.

- `--MinEditSites 2`: Set minimum number of editing sites per gene. Only genes containing at least 2 editing sites will be reported.

- `--EditedMinCoverage 5`: Set minimum coverage requirement for the edited group. Sites in the edited group must have coverage of 5 or above.

- `--ControlMinCoverage 5`: Set minimum coverage requirement for the control group. Sites in the control group must have coverage of 5 or above.

- `--cpu 8`: Set the number of CPU cores for parallel processing.

- `--outfile "$OUTFILE"`: Specify the output result file path.

- `--verbose`: Enable verbose output mode to display detailed information during processing.

### 8. Merge All Results

Merge all results in batches using a two-step merging strategy: first merge results from individual wild-type and all mutant types, then merge all wild-type final results.

#### Step 1: Merge Individual Wild-type with All Mutant Results

```bash
perl summarize_sites.pl --repOnly "${MUT_BED_FILES[@]}" > "$WT_OUTPUT_FILE"
```

#### Step 2: Merge All Wild-type Final Results

```bash
perl summarize_sites.pl --minRep all*60% "${WT_OUTPUT_FILES[@]}" > "$FINAL_OUTPUT_FILE"
```

#### Parameter Description

##### First Step Merging Parameters

- `--repOnly`: Only report editing sites that appear repeatedly in specified files, used to filter reliable editing sites that appear in multiple mutant samples.

- `"${MUT_BED_FILES[@]}"`: Array of BED format result files from all mutant samples.

- `> "$WT_OUTPUT_FILE"`: Output merged results from individual wild-type vs mutant comparison.

##### Second Step Merging Parameters

- `--minRep all*60%`: Set minimum repetition rate threshold to 60%. Only editing sites that appear in at least 60% of wild-type samples will be included in the final results, ensuring reliability and consistency.

- `"${WT_OUTPUT_FILES[@]}"`: Array of all wild-type sample processed output files.

- `> "$FINAL_OUTPUT_FILE"`: Final merged output result file.

#### Processing Workflow Description

1. **First Phase**: For each wild-type sample, compare it with all mutant samples to filter editing sites that repeatedly appear in multiple mutants.

2. **Second Phase**: Aggregate results from all wild-type samples, retaining only high-confidence editing sites that can be detected in the majority (≥60%) of wild-type samples.

This hierarchical merging strategy ensures high quality and biological significance of the final results.

At this point, we have obtained the methylation count data.











