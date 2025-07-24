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
Note: It is strongly recommended to use the specified software versions to avoid compatibility issues.

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

## 













