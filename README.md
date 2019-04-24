## ChIPseqSpikeInFree

A Spike-in Free ChIP-Seq Normalization Approach for Detecting Global Changes in Histone Modifications

## Background

Global histone modification changes can be detected by using exogenous reference spike-in controls. However, much of the ChIP-seq data available in public repositories do not contain spike-in controls. To make quantitative comparisons between data with and without spike-in, researchers must regenerate whole datasets with spike-in ChIP-seq protocols. However, this is not always a feasible solution. A basic scaling factor calculation for these scenarios remains a challenge.

We introduce `ChIPseqSpikeInFree`, a novel ChIP-seq normalization method to effectively determine scaling factors for samples across various conditions and treatments, which does not rely on exogenous spike-in chromatin or peak detection to reveal global changes in histone modification occupancy. This method is capable of revealing the same magnitude of global changes as the spike-in method.

## Prerequisites

`ChIPseqSpikeInFree` depends on `Rsamtools`, `GenomicRanges`, and `GenomicAlignments` to count reads from bam files.

To install these packages, start `R` (version "3.4") and enter:
```R
> source("https://bioconductor.org/biocLite.R")
> biocLite("Rsamtools")
> biocLite("GenomicRanges")
> biocLite("GenomicAlignments")
```

If you use `R` (version "3.5"), enter:
```R
> if (!requireNamespace("BiocManager", quietly = TRUE))
>     install.packages("BiocManager")
> BiocManager::install("Rsamtools")
> BiocManager::install("GenomicRanges")
> BiocManager::install("GenomicAlignments")
```

## Installation

Using `R`, enter:

```R
# Option 1. Install this package from CRAN
> install.packages("ChIPseqSpikeInFree")

# Option 2. Install this package from GitHub
> install.packages("devtools")
> library(devtools)
> install_github("hongjianjin/ChIPseqSpikeInFree")
```

### Usage

A simple workflow in `R` environment.

##### 0. Load package

```R
> library("ChIPseqSpikeInFree")
```

##### 1. Generate `sample_meta.txt` (**tab-delimited txt file with header line**) as follows

Save as `/your/path/sample_meta.txt`

| ID           | ANTIBODY | GROUP |
| ------------ | -------- | ----- |
| ChIPseq1.bam | H3K27me3 | WT    |
| ChIPseq2.bam | H3K27me3 | K27M  |

```R
> metaFile <- "/your/path/sample_meta.txt"
```

##### 2. Assign bam file names to a vector

```R
> bams <- c("ChIPseq1.bam","ChIPseq2.bam")
```

##### 3. Run `ChIPseqSpikeInFree` pipeline (when your bam files correspond to the human reference hg19)

```R
> ChIPseqSpikeInFree(bamFiles = bams, chromFile = "hg19", metaFile = metaFile, prefix = "test")
```

### Input

In the simple usage scenario, the user should have ChIP-seq bam files ready. Sample information can be specified in a metadata file (`metaFile`) and the user should choose a correct reference genome corresponding to the bams.

##### 1. bamFiles: a vector of bam filenames

User should follow ChIP-seq guidelines suggested by `ENCODE consortium(Landt, et al., 2012)` and check the data quality first. We highly recommend the following:

1) **Remove low-quality or non-unique reads**
2) **Remove spike-in reads** from your bam files before you run `ChIPseqSpikeInFree` normalization
3) Your bam files must **contain a header section** and an alignment section.

##### 2. chromFile: chromosome sizes of reference genome.
`hg19`, `mm9`, `mm10`, `hg19` are included in the package.

For other genomes, you can either:

* Use `fetchChromSizes` to get it from UCSC, but not all genomes are available.

    `http://hgdownload.soe.ucsc.edu/goldenPath/${DB}/bigZips/${DB}.chrom.sizes`  
    Replace `${DB}` with reference genome

* Generate directly from `fasta` file (Linux)
    ```bash
    $ samtools faidx genome.fa
    $ cut -f1,2 genome.fa.fai > genome.chrom.sizes
    ```

##### 3. metaFile:
A tab-delimited text file **having three columns with a header line: `ID`, `ANTIBODY`, and `GROUP`**.

* `ID` is the bam file name of ChIP-seq sample that will be included for analysis.
* `ANTIBODY` represents antibody used for ChIP.
* `GROUP` describes the biological treatment or condition of this sample.

### Output

After you successfully run following `ChIPseqSpikeInFree` pipeline:
```R
> ChIPseqSpikeInFree(bamFiles = bams, chromFile = "hg19", metaFile = metaFile, prefix = "test")
```

Output will include: (in case that you set `prefix ="test"`)
1. `test_SF.txt` - text result
    * tab-delimited text format, a table of calculated scaling factors by pipeline
2. `test_boxplot.pdf` - graphical result
    * view of scaling factors as boxplot based on `test_SF.txt`
3. `test_rawCounts.txt` - intermediate file
    * tab-delimited text format, a table of raw read counts for each 1kb bin across genome
4. `test_parsedMatrix.txt` - intermediate file
    * tab-delimited text format, a table of proportion of reads below given cutoffs (CPMW)
5. `test_distribution.pdf` - intermediate file
    * view of proportion of reads below the given CPMW based on `test_parsedMatrix.txt`

## Notes

This repository contains the following:
- source code
- documentation
- chromFile of human and mouse reference genome (`hg19`, `mm9`, `mm10`, and `hg38`)
- an example of `sample_meta.txt`

Users will still need to source and provide the following:
- bam files
- `sample_meta.txt`
- chromFile except `hg19`, `mm9`, `mm10`, and `hg38`

## Maintainers

* [Hongjian Jin]

[Hongjian Jin]: https://github.com/hongjianjin
