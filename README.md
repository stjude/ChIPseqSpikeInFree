## ChIPseqSpikeFree, 
A Spike-in Free ChIP-Seq Normalization Approach for Detecting Global Changes in Histone Modifications

## Background

The detection of global histone modification changes can be addressed using exogenous reference Spike-in controls. However, many thousands of ChIP-seq data available in public depositories nowadays were done without including Spike-in procedure. In order to do quantitative comparisons between these data, researchers have to regenerate whole data set using spikein ChIP-seq protocols – this is an infeasible solution sometime. A basic scaling factor calculation for these scenarios remains a problem with surprisingly few solutions presented so far. We pesent ChIPseqSpikeFree , a simple ChIP-seq normalization method to effectively determine scaling factors for samples across different conditions or treatments, which doesn't rely on exogenous spike-in chromatin or peak detection to reveal global changes in histone modification occupancy. It can reveal same magnitude of global changes compared to spike-In method.

## Prerequisites

ChIPseqSpikeFree depends on Rsamtools,GenomicRanges and GenomicAlignments to count reads from bam file.
To install these packages, start R (version "3.4") and enter:
```
>source("https://bioconductor.org/biocLite.R")
>biocLite("Rsamtools")
>biocLite("GenomicRanges")
>biocLite("GenomicAlignments")
```
If you use R (version "3.5") and enter:
```
>if (!requireNamespace("BiocManager", quietly = TRUE))
>    install.packages("BiocManager")
>BiocManager::install("Rsamtools")
>BiocManager::install("GenomicRanges")
>BiocManager::install("GenomicAlignments")
```

## Installation

If you use R, enter
```
#Option 1. intall this package from CRAN
>install.packages("ChIPseqSpikeFree")

#Option 2. intall this package from GitHub
>install.packages("devtools")
>library(devtools)
>install_github("hongjianjin/ChIPseqSpikeFree")
```

### Usage

A simple workflow.
```
##0. load package
>library("ChIPseqSpikeFree")
```
##1. generate a sample_meta.txt (tab-delimited txt file) as follows
#/your/path/sample_meta.txt
# ID ANTIBODY GROUP
# ChIPseq1.bam H3K27me3 WT
# ChIPseq2.bam H3K27me3 K27M

```
>metaFile <- "/your/path/sample_meta.txt"
>meta <- ReadMeta(metaFile)
>head(meta)
ID ANTIBODY GROUP
ChIPseq1.bam H3K27me3 WT
ChIPseq2.bam H3K27me3 K27M

##2. assign bam file names to a vector
>bams <- c("ChIPseq1.bam","ChIPseq2.bam")
##3. run ChIPseqSpikeFree pipeline (when your bam files correspond to the human reference hg19) 
>ChIPseqSpikeFree(bamFiles=bams, chromFile="hg19",metaFile=metaFile,prefix="test")
```

### Input

In the simple usage scenario, the user should have ChIP-seq Bam files ready and sample information can be specified in a metadata file (metaFile) and should choose a correct reference genome corresponding to bams. 

##### 1.bamFiles: a vector of bam filenames

User should follow ChIP-seq guidelines suggested by EMCODE consortium(Landt, et al., 2012) and check the data quality. We recommend to remove low-quality or non-unique reads from your bam files before you run ChIPseqSpikeFree normalization.

##### 2.chromFile:chromosome sizes of reference genome. 
"hg19", "mm9","mm10","hg19" are included in the package.
For other genomes, you could 
- 2.1 use fetchChromSizes to get it from UCSC, but not all genomes are available. (replace following ${DB} with reference genome)
"http://hgdownload.soe.ucsc.edu/goldenPath/${DB}/bigZips/${DB}.chrom.sizes"
- 2.2 generate directly from fasta file (Linux)
$samtools faidx genome.fa
$cut -f1,2 genome.fa.fai > genome.chrom.sizes

##### 3.metaFile: 
A tab-delimited text file having three columns: ID, ANTIBODY and GROUP. Where ID is the bam file name of ChIP-seq sample that will be included for analysis;  ANTIBODY represents antibody used for ChIP and GROUP describes the biological treatment or condition of this sample. 


### Output

After you successfully run following ChIPseqSpikeFree pipeline 
```
>ChIPseqSpikeFree(bamFiles=bams, chromFile="hg19",metaFile=metaFile,prefix="test")
```
Output will include: (in case that you set prefix ="test")
##### 1. test_SF.txt (text result)
- tab-delimited text formart, a table of caculated scaling factors by pipeline
##### 2. test_boxplot.pdf (graphical result)
- view of scaling factors as boxplot based on test_SF.txt
##### 3. test_rawCounts.txt (intermediate file)
- tab-delimited text formart, a table of raw read count for each 1kb bin across genome
##### 4. test_parsedMatrix.txt (intermediate file)
- tab-delimited text formart, a table of proportion of reads below given cutoffs (CPMW)
##### 5. test_distribution.pdf (intermediate file)
- view of proportion of reads below the given CPMW based on test_parsedMatrix.txt


## Notes

1.What **IS** included in this repository 
- source code
- documentation
- chromFile of human and mouse reference genome (hg19, mm9, mm10 and hg38)
- an example of sample_meta.txt

2.What is **NOT** included and users will have to provide by themselves:
- bam files
- sample_meta.txt
- chromFile except hg19, mm9, mm10 and hg38

