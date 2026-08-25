---
name: chipseqspikeinfree
description: Run the ChIPseqSpikeInFree command-line workflow for spike-in-free ChIP-seq normalization from BAM files, metadata, and genome chromosome sizes; also support replotting from an existing parsed matrix.
usage: Use when a user wants an AI agent to validate ChIPseqSpikeInFree inputs, choose standard or replotting mode, run the CLI wrapper, and report scaling-factor outputs and warnings.
version: 1.1.0
validated_against:
  date: 2026-08-25
  package: ChIPseqSpikeInFree
  package_version: 1.2.6
  cli_readme: ChIPseqSpikeInFree_cli_README.md
  cli_script: inst/scripts/run_ChIPseqSpikeInFree.R
tags:
  - skill
  - category:pipeline
  - domain:genomics
  - domain:chip-seq
  - language:r
  - interface:cli
---

# ChIPseqSpikeInFree

## Skill Overview

### [ChIPseqSpikeInFree](https://github.com/stjude/ChIPseqSpikeInFree)

Maintainer: Hongjian Jin, St. Jude Children's Research Hospital.

### Purpose

This skill enables an AI agent to run the ChIPseqSpikeInFree command-line workflow for spike-in-free ChIP-seq normalization. The agent validates user inputs, selects standard or replotting mode, runs the CLI wrapper, and summarizes scaling-factor outputs, QC status, and any warnings or failures.

## Usage For The AI Agent

Use this skill when the user asks to run ChIPseqSpikeInFree from the command line, calculate spike-in-free scaling factors from BAM files, regenerate plots from an existing parsed matrix, prepare CLI inputs, or interpret CLI outputs.

The agent should:

1. Identify whether the user wants standard mode or replotting mode.
2. Collect and validate the required inputs for that mode.
3. Locate the ChIPseqSpikeInFree repository or use the current repository checkout when already inside it.
4. Confirm that R, the `ChIPseqSpikeInFree` package, `optparse`, and required Bioconductor dependencies are available.
5. Build a command using absolute paths for user-supplied files and output prefixes.
6. Run the workflow in the requested output location without overwriting existing results unless the user explicitly asks for replacement.
7. Report the command, output prefix, selected mode, key parameters, generated files, and warnings or failures.

Prefer the CLI wrapper documented in `ChIPseqSpikeInFree_cli_README.md`:

```bash
Rscript inst/scripts/run_ChIPseqSpikeInFree.R \
  -L /absolute/path/bam.lst \
  -G hg19 \
  -m /absolute/path/sample_meta.txt \
  -o /absolute/path/output/test
```

## Information The User Needs To Provide

### Required For Standard Mode

Standard mode counts reads, parses count distributions, calculates scaling factors, and generates plots.

Required inputs:

- `bamList`: plain text file with one BAM path per line.
- `genome`: one of `hg19`, `hg38`, `mm9`, `mm10`, or an absolute path to a chromosome-size file.
- `metaData`: tab-delimited metadata file with `ID`, `ANTIBODY`, and `GROUP` columns.
- `prefix`: output prefix, preferably an absolute path such as `/path/to/results/test`.

Standard-mode command template:

```bash
Rscript inst/scripts/run_ChIPseqSpikeInFree.R \
  -L <bamList> \
  -G <genome-or-chrom-size-file> \
  -m <metadata-file> \
  -o <output-prefix>
```

### Required For Replotting Mode

Replotting mode reuses an existing parsed matrix and regenerates scaling-factor and plot outputs.

Required inputs:

- `metaData`: tab-delimited metadata file with `ID`, `ANTIBODY`, and `GROUP` columns.
- `prefix`: the original output prefix whose `<prefix>_parsedMatrix.txt` already exists.
- `replot`: must be `TRUE`.

Replotting command template:

```bash
Rscript inst/scripts/run_ChIPseqSpikeInFree.R \
  -m <metadata-file> \
  -o <output-prefix> \
  -r TRUE
```

Read-counting options such as `-L`, `-G`, `-b`, `-E`, and `-n` are ignored in replotting mode.

### Optional Inputs

Option | Accepted Values | Default Behavior
--- | --- | ---
`binSize` | Numeric bp value from 200 to 10000 | `1000`
`paired` | `TRUE` or `FALSE` | `FALSE`, meaning single-end counting
`minFirstTurn` | Numeric fraction from about 0.20 to 0.80, or omit | automatic first-turn selection
`maxLastTurn` | Numeric fraction, commonly 0.90 to 0.99 | `0.99`
`cutoffQC` | Numeric QC cutoff, commonly 0.95 to 1.2 | `1.2`
`xmax` | Numeric x-axis maximum for distribution plots | automatic plot limit
`ncores` | Positive integer | `2`
`checkout_directory` | Existing or new repository checkout path | current repository if already present
`output_directory` | Existing or new results directory | infer from `prefix` when safe

## Input File Rules

### BAM List

The BAM list must contain one BAM path per line. Confirm that every BAM file exists and is readable before running standard mode.

```text
/absolute/path/sample1.bam
/absolute/path/sample2.bam
/absolute/path/sample3.bam
```

### Metadata

The metadata file must be tab-delimited and include a header with at least:

```text
ID    ANTIBODY    GROUP
sample1.bam    H3K27me3    WT
sample2.bam    H3K27me3    Treatment
```

Rules:

- `ID` values must be unique.
- `ID` values must match BAM basenames or parsed count matrix column names.
- `COLOR` is optional and used for plots when present.
- Samples that will be compared should be processed together in the same run and antibody context.

### Genome Or Chromosome Sizes

Use a built-in genome alias when possible:

```text
hg19
hg38
mm9
mm10
```

For other genomes, provide a readable chromosome-size file with chromosome name and size columns.

## Validation And Safety Rules

Before running the workflow, the agent should:

- Confirm that the repository checkout contains `inst/scripts/run_ChIPseqSpikeInFree.R`.
- Confirm that required paths exist and are readable.
- Confirm that the output directory exists or can be created safely.
- Check whether `<prefix>_rawCounts.txt`, `<prefix>_parsedMatrix.txt`, `<prefix>_SF.txt`, `<prefix>_distribution.pdf`, or `<prefix>_boxplot.pdf` already exist.
- Ask before deleting, replacing, or intentionally overwriting existing outputs.
- Validate that `binSize` is between 200 and 10000 when supplied.
- Validate that paired-end mode is requested only when the BAM files contain paired-end alignments.
- In replotting mode, confirm that `<prefix>_parsedMatrix.txt` exists before running.
- Surface missing, ambiguous, or contradictory information before execution when it cannot be resolved safely.

Use absolute paths in commands whenever possible. Clearly distinguish input files from generated output files.

## Output Files

For `prefix = test`, standard mode writes:

```text
test_rawCounts.txt
test_parsedMatrix.txt
test_SF.txt
test_distribution.pdf
test_boxplot.pdf
```

Replotting mode writes scaling-factor and plot outputs with a `_replot` suffix, such as:

```text
test_replot_SF.txt
test_replot_distribution.pdf
test_replot_boxplot.pdf
```

The CLI may reuse existing `<prefix>_rawCounts.txt` or `<prefix>_parsedMatrix.txt` to avoid repeating expensive steps. Tell the user when existing intermediates are reused.

## Interpretation Guidance

- `QC = pass` means the sample is eligible for scaling-factor interpretation.
- Failed-QC samples should not be interpreted as valid scaling-factor estimates.
- `SF` is calculated within each antibody group.
- Larger `SF` indicates lower global histone mark signal relative to the run-specific reference sample.
- Scaling factors should not be compared across independent runs unless all relevant samples were processed together.

## Recommended Agent Prompt

```text
Use $chipseqspikeinfree.
Locate or clone ChIPseqSpikeInFree in <checkout-directory>.
My BAM list is <bam-list-file>.
My metadata file is <metadata-file>.
Use genome <hg19|hg38|mm9|mm10|chrom-size-file>.
Write outputs with prefix <output-prefix>.
Use standard mode.
Optional settings: binSize=<value>, paired=<TRUE|FALSE>, maxLastTurn=<value>, cutoffQC=<value>, ncores=<value>.
```

For replotting:

```text
Use $chipseqspikeinfree.
Locate ChIPseqSpikeInFree in <checkout-directory>.
My metadata file is <metadata-file>.
The existing output prefix is <output-prefix> and <output-prefix>_parsedMatrix.txt already exists.
Use replotting mode.
Optional settings: minFirstTurn=<value>, maxLastTurn=<value>, cutoffQC=<value>, xmax=<value>.
```

## Expected Agent Response

After execution, summarize:

- whether the workflow completed successfully
- the exact command that was run
- the selected mode and key parameters
- the output prefix and result directory
- whether existing raw-count or parsed-matrix intermediates were reused
- which expected files were produced
- any warnings, missing outputs, or failed samples requiring review
- the next useful action, such as inspecting `*_SF.txt` or opening the generated PDFs

## Maintenance Notes

Keep this skill platform-neutral. Do not add Codex-specific, Claude-specific, Gemini-specific, or local-machine-only instructions to the core skill. If the CLI behavior changes, update this skill and `ChIPseqSpikeInFree_cli_README.md` together.

Validate package tests from the repository root with:

```bash
Rscript -e "testthat::test_dir('tests/testthat')"
```

Success criterion: `FAIL 0`. Review warnings, but known locale warnings from constrained environments do not necessarily indicate workflow failure.
