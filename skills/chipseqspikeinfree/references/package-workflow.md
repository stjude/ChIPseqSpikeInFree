# ChIPseqSpikeInFree Package Workflow

## Repository Shape

- Main source: `R/ChIPseqSpikeInFree.R`
- Package metadata: `DESCRIPTION`
- Exports/imports: `NAMESPACE`
- User docs: `README.md`, `man/*.Rd`, `docs/`
- Test suite: `tests/testthat/`
- Current package version: `1.2.6`

## Inputs

Metadata files must be tab-delimited with a header and at least:

- `ID`: BAM basename matching count matrix/sample names
- `ANTIBODY`: ChIP antibody or mark
- `GROUP`: biological condition/treatment
- `COLOR`: optional plotting color

BAM files should be quality-filtered, duplicate-filtered if appropriate, indexed when needed, and aligned to the same genome as `chromFile`.

Supported bundled chromosome-size aliases are `hg19`, `hg38`, `mm9`, and `mm10`. Custom chromosome-size files are tab-delimited `chrom size` files.

## Main Functions

- `GenerateBins(chromFile, binSize, overlap, withChr, minChrSize)`: create genome bins. `binSize` should be 200-10000 bp; `overlap` must be >= 0 and smaller than `binSize`.
- `CountRawReads(bamFiles, chromFile, prefix, singleEnd, binSize)`: count reads per bin and write `<prefix>_rawCounts.txt`.
- `ReadMeta(metaFile)`: read and validate metadata; duplicate IDs are errors.
- `ParseReadCounts(data, metaFile, by, prefix, binSize, ncores)`: write `<prefix>_parsedMatrix.txt`.
- `CalculateSF(data, metaFile, minFirstTurn, maxLastTurn, cutoff_QC)`: return metadata with QC, turns, slopes, and `SF`; failed-QC samples have `SF = NA`.
- `PlotDistr(data, SF, prefix, xlimMaxCPMW)`: write `<prefix>_distribution.pdf` and `<prefix>_SF.txt`.
- `BoxplotSF(input, prefix)`: write `<prefix>_boxplot.pdf`.
- `ChIPseqSpikeInFree(...)`: wrapper; invisibly returns the scaling-factor metadata table.

## Outputs

For `prefix = "test"`, expected outputs are:

- `test_rawCounts.txt`
- `test_parsedMatrix.txt`
- `test_SF.txt`
- `test_distribution.pdf`
- `test_boxplot.pdf`

Interpret larger scaling factors as lower global histone mark signal relative to the selected reference within the same antibody. Do not compare scaling factors across independently processed batches.

## Testing

Preferred command from repo root:

```powershell
powershell -ExecutionPolicy Bypass -File skills\chipseqspikeinfree\scripts\run_tests.ps1
```

Equivalent direct command on this Windows environment:

```powershell
& 'C:\Users\hjin\AppData\Local\Programs\R\R-4.3.3\bin\x64\Rscript.exe' -e "testthat::test_dir('tests/testthat')"
```

Expected current focused-suite result: `FAIL 0`, `WARN 0`, `SKIP 0`, `PASS 17`.

## Release Notes For Agents

- Keep README version snippets in sync with `DESCRIPTION`.
- Keep the README PDF link pointing to an existing file in `docs/`; currently only `docs/ChIPseqSpikeInFree_1.2.4.pdf` exists.
- This repo may be on Windows/PowerShell and `Rscript` may not be on PATH. Use the fallback path above when needed.
- Avoid changing generated `.Rd` files unless regenerating documentation is part of the user request.