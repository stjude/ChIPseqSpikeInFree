# ChIPseqSpikeInFree CLI

`inst/scripts/run_ChIPseqSpikeInFree.R` is a command-line wrapper for running the ChIPseqSpikeInFree spike-in-free ChIP-seq normalization workflow from BAM files, metadata, and a genome chromosome-size source.

The CLI supports two modes:

- Standard mode: count reads, parse count distributions, calculate scaling factors, and generate plots.
- Replotting mode: reuse an existing parsed matrix to recalculate scaling factors and regenerate plots.

## Requirements

- R 3.6 or newer.
- Installed R package: `ChIPseqSpikeInFree`.
- Installed R package: `optparse`.
- ChIPseqSpikeInFree runtime dependencies, including Bioconductor packages used for BAM reading and genomic interval counting.
- Input BAM files aligned to the same genome build.

The script exits early if the detected R version is older than 3.6.

## Input Files

### BAM List

The BAM list is a plain text file with one BAM path per line:

```text
/path/to/sample1.bam
/path/to/sample2.bam
/path/to/sample3.bam
```

Use the same BAM files represented in the metadata file. The CLI uses BAM basenames as count matrix column names, so metadata `ID` values should match the BAM basenames.

### Metadata

The metadata file must be tab-delimited and include at least these columns:

```text
ID    ANTIBODY    GROUP
sample1.bam    H3K27me3    WT
sample2.bam    H3K27me3    Treatment
```

Optional column:

- `COLOR`: used for plotting when present.

Rules:

- `ID` values must be unique.
- `ID` values must match BAM basenames or parsed count matrix column names.
- Scaling factors are only comparable within the same jointly processed run and antibody context.

### Genome Or Chromosome Sizes

Use one of the built-in genome aliases:

```text
hg19
hg38
mm9
mm10
```

or provide a custom chromosome-size file.

## Standard Mode

Standard mode requires:

- `-L`, `--bamList`
- `-G`, `--genome`
- `-m`, `--metaData`
- `-o`, `--prefix`

Example:

```bash
Rscript inst/scripts/run_ChIPseqSpikeInFree.R \
  -L bam.lst \
  -G hg19 \
  -m sample_meta.txt \
  -o test
```

Example for a mark such as H3K9me3, where a lower QC cutoff and lower `maxLastTurn` may be useful:

```bash
Rscript inst/scripts/run_ChIPseqSpikeInFree.R \
  -L bam.lst \
  -G hg19 \
  -m sample_meta.txt \
  -o H3K9me3 \
  -l 0.97 \
  -c 1
```

## Replotting Mode

Replotting mode uses an existing parsed matrix:

```text
<prefix>_parsedMatrix.txt
```

It recalculates scaling factors and writes outputs with a `_replot` suffix. Read-counting options such as `-L`, `-G`, `-b`, `-E`, and `-n` are ignored in this mode.

Example:

```bash
Rscript inst/scripts/run_ChIPseqSpikeInFree.R \
  -m sample_meta.txt \
  -o test \
  -r TRUE
```

Example with a custom x-axis maximum for the distribution plot:

```bash
Rscript inst/scripts/run_ChIPseqSpikeInFree.R \
  -m sample_meta.txt \
  -o test \
  -x 10 \
  -r TRUE
```

## Options

Option | Meaning | Default
--- | --- | ---
`-L`, `--bamList` | Text file with one BAM path per line. Required in standard mode. | `NA`
`-G`, `--genome` | Genome alias or chromosome-size file. Required in standard mode. | `NA`
`-m`, `--metaData` | Metadata file. Required in both modes. | `NA`
`-o`, `--prefix` | Output prefix. Required in both modes. | `NA`
`-b`, `--binSize` | Bin size in bp for read counting. | `1000`
`-r`, `--replot` | Set to `TRUE` to run replotting mode. | `FALSE`
`-E`, `--paired` | Set to `TRUE` for paired-end counting. The default is single-end counting. | `FALSE`
`-f`, `--minFirstTurn` | Lower fraction for first turn-point estimation. Uses automatic selection when omitted. | `NA`
`-l`, `--maxLastTurn` | Upper fraction used to estimate the last turn point. | `0.99`
`-c`, `--cutoffQC` | QC cutoff. Smaller values are looser. | `1.2`
`-x`, `--xmax` | Maximum x-axis value for distribution plots. | `NULL`
`-n`, `--ncores` | Number of worker cores for parsing raw counts. | `2`

## Expected Outputs

For `prefix = test`, standard mode writes:

```text
test_rawCounts.txt
test_parsedMatrix.txt
test_SF.txt
test_distribution.pdf
test_boxplot.pdf
```

Replotting mode writes plot and scaling-factor outputs with the prefix `test_replot`.

## Resume Behavior

The CLI avoids repeating expensive steps when intermediate files already exist:

- If `<prefix>_rawCounts.txt` exists, the script loads it instead of recounting.
- If `<prefix>_parsedMatrix.txt` exists, the script loads it instead of reparsing.

Delete those files before rerunning if you need to force recounting or reparsing.

## Interpretation

- `QC = pass` means the sample is eligible for scaling-factor interpretation.
- Failed-QC samples should not be interpreted as valid scaling-factor estimates.
- `SF` is calculated within each antibody group.
- Larger `SF` indicates lower global histone mark signal relative to the run-specific reference sample.
- Do not compare scaling factors from independent runs unless all relevant samples were processed together.

## Implementation Notes

The CLI keeps orchestration logic in `ChIPseqSpikeInFree2()` and delegates the core workflow steps to the exported package functions where possible. This keeps command-line behavior aligned with the maintained implementations in `R/ChIPseqSpikeInFree.R`.

Current behavior:

- `--paired TRUE` runs paired-end counting; the default is single-end counting.
- `--minFirstTurn` is passed to `CalculateSF()` in both standard and replotting modes. When omitted, the CLI uses automatic first-turn selection.
- Out-of-range `binSize` values below 200 bp or above 10000 bp are rejected before the workflow runs.
- The wrapper invisibly returns the scaling-factor table from `ChIPseqSpikeInFree2()`.

For new development, keep CLI behavior aligned with the exported package functions and add tests for command-line behavior before changing option semantics.

## Validation

Run the package test suite from the repository root:

```bash
Rscript -e "testthat::test_dir('tests/testthat')"
```

Success criterion:

```text
FAIL 0
```


