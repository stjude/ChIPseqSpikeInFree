# ChIPseqSpikeInFree Skill

An AI-agent skill for running the [ChIPseqSpikeInFree](https://github.com/stjude/ChIPseqSpikeInFree) spike-in-free ChIP-seq normalization workflow.

This skill is designed for skill-aware AI agents. It teaches the agent how to:

- validate BAM lists, metadata files, genome aliases, and chromosome-size files
- choose standard mode or replotting mode
- run the `inst/scripts/run_ChIPseqSpikeInFree.R` CLI wrapper
- reuse existing raw-count or parsed-matrix intermediates safely
- inspect generated scaling-factor tables and plots
- summarize QC status, scaling factors, warnings, and output locations

## What The User Needs To Provide

The user does not need to understand the ChIPseqSpikeInFree repository internals. They only need to provide the workflow inputs for the mode they want to run.

Standard mode requires:

- a BAM list file with one BAM path per line
- a metadata file with `ID`, `ANTIBODY`, and `GROUP` columns
- a genome alias (`hg19`, `hg38`, `mm9`, or `mm10`) or a chromosome-size file
- an output prefix

Replotting mode requires:

- a metadata file with `ID`, `ANTIBODY`, and `GROUP` columns
- an existing output prefix whose `<prefix>_parsedMatrix.txt` file already exists

Optional settings include:

- `binSize`, default `1000`, accepted range `200` to `10000`
- `paired`, default `FALSE`, meaning single-end counting
- `minFirstTurn`, default automatic first-turn selection
- `maxLastTurn`, default `0.99`
- `cutoffQC`, default `1.2`
- `xmax`, default automatic plot limit
- `ncores`, default `2`

## Recommended Agent Prompt

For standard mode:

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

For replotting mode:

```text
Use $chipseqspikeinfree.
Locate ChIPseqSpikeInFree in <checkout-directory>.
My metadata file is <metadata-file>.
The existing output prefix is <output-prefix> and <output-prefix>_parsedMatrix.txt already exists.
Use replotting mode.
Optional settings: minFirstTurn=<value>, maxLastTurn=<value>, cutoffQC=<value>, xmax=<value>.
```

## What The Skill Does

When invoked, the skill instructs the agent to:

1. locate the ChIPseqSpikeInFree repository or use the current checkout
2. confirm the CLI wrapper exists at `inst/scripts/run_ChIPseqSpikeInFree.R`
3. validate required paths and metadata before execution
4. select standard mode or replotting mode from the user's request
5. check for existing output files before running
6. build a command with absolute paths whenever possible
7. run the CLI workflow
8. confirm expected outputs were produced
9. report QC status, scaling factors, output files, and warnings

The skill should preserve existing user files unless the user explicitly requests replacement.

## Direct Command-Line Use

You can also run the ChIPseqSpikeInFree CLI directly without going through an agent.

Standard mode example:

```bash
Rscript inst/scripts/run_ChIPseqSpikeInFree.R \
  -L /path/to/bam.lst \
  -G hg19 \
  -m /path/to/sample_meta.txt \
  -o /path/to/results/test
```

Standard mode with broad-mark settings:

```bash
Rscript inst/scripts/run_ChIPseqSpikeInFree.R \
  -L /path/to/bam.lst \
  -G hg19 \
  -m /path/to/sample_meta.txt \
  -o /path/to/results/H3K9me3 \
  -l 0.97 \
  -c 1
```

Paired-end standard mode:

```bash
Rscript inst/scripts/run_ChIPseqSpikeInFree.R \
  -L /path/to/bam.lst \
  -G hg19 \
  -m /path/to/sample_meta.txt \
  -o /path/to/results/test \
  -E TRUE
```

Replotting mode example:

```bash
Rscript inst/scripts/run_ChIPseqSpikeInFree.R \
  -m /path/to/sample_meta.txt \
  -o /path/to/results/test \
  -r TRUE
```

Replotting mode with a custom x-axis maximum:

```bash
Rscript inst/scripts/run_ChIPseqSpikeInFree.R \
  -m /path/to/sample_meta.txt \
  -o /path/to/results/test \
  -x 10 \
  -r TRUE
```

## Expected Outputs

For output prefix `test`, standard mode writes:

- `test_rawCounts.txt`
- `test_parsedMatrix.txt`
- `test_SF.txt`
- `test_distribution.pdf`
- `test_boxplot.pdf`

Replotting mode writes scaling-factor and plot outputs with a `_replot` suffix, such as:

- `test_replot_SF.txt`
- `test_replot_distribution.pdf`
- `test_replot_boxplot.pdf`

The CLI may reuse existing `test_rawCounts.txt` or `test_parsedMatrix.txt` to avoid repeating expensive steps. Delete those intermediates first only when recounting or reparsing is intended.

## Interpreting Results

The most important output is `<prefix>_SF.txt`.

Key columns:

- `QC`: sample quality-control status
- `SF`: scaling factor for passing samples
- `TURNS`: turn-point coordinates used for slope-based scaling-factor calculation

Interpretation rules:

- `QC = pass` means the sample is eligible for scaling-factor interpretation.
- Failed-QC samples should not be treated as valid scaling-factor estimates.
- Scaling factors are calculated within each antibody group.
- Larger `SF` indicates lower global histone mark signal relative to the run-specific reference sample.
- Scaling factors should not be compared across independent runs unless all relevant samples were processed together.

## Behavior Verified In Testing

This skill has been aligned with:

- `ChIPseqSpikeInFree` package version `1.2.6`
- `inst/scripts/run_ChIPseqSpikeInFree.R`
- `ChIPseqSpikeInFree_cli_README.md`
- the focused `testthat` suite in `tests/testthat`

Validation command:

```bash
Rscript -e "testthat::test_dir('tests/testthat')"
```

Expected success criterion:

```text
FAIL 0
```

Warnings should be reviewed. Known locale warnings from constrained environments do not necessarily indicate package or workflow failure.

## Troubleshooting

- If `Rscript` is missing, install R or provide the full path to `Rscript`.
- If `optparse` is missing, install the R package before running the CLI.
- If Bioconductor dependencies are missing, install the package dependencies for `ChIPseqSpikeInFree`.
- If a BAM path fails validation, check the BAM list and use absolute paths.
- If metadata IDs do not match count matrix columns, make sure `ID` values match BAM basenames.
- If replotting mode fails, confirm that `<prefix>_parsedMatrix.txt` exists.
- If expected outputs already exist, decide whether to reuse them or remove them before rerunning.
- If no samples pass QC, inspect enrichment quality, antibody grouping, `maxLastTurn`, and `cutoffQC`.

## Files In This Skill

- [`SKILL.md`](SKILL.md): core platform-neutral skill instructions for AI agents
- [`README.md`](README.md): human-readable overview and quick-start guide for the skill

Related package files:

- [`../../ChIPseqSpikeInFree_cli_README.md`](../../ChIPseqSpikeInFree_cli_README.md): detailed CLI documentation
- [`../../inst/scripts/run_ChIPseqSpikeInFree.R`](../../inst/scripts/run_ChIPseqSpikeInFree.R): command-line wrapper
- [`../../README.md`](../../README.md): main package README
