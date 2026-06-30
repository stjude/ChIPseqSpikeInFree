test_that("GenerateBins validates bin size and overlap", {
  chrom_file <- tempfile(fileext = ".chrom.sizes")
  writeLines("chr1\t1000", chrom_file)

  expect_error(GenerateBins(chrom_file, binSize = 199), "Recommended binSize")
  expect_error(GenerateBins(chrom_file, binSize = 10001), "Recommended binSize")
  expect_error(GenerateBins(chrom_file, binSize = 200, overlap = 200), "overlap")
  expect_error(GenerateBins(chrom_file, binSize = 200, overlap = -1), "overlap")
})

test_that("GenerateBins defaults minChrSize before reporting excluded contigs", {
  chrom_file <- tempfile(fileext = ".chrom.sizes")
  writeLines(c("chr1\t500", "chrUn\t1000"), chrom_file)

  bins <- expect_warning(
    GenerateBins(chrom_file, binSize = 200, overlap = 0, minChrSize = NULL),
    NA
  )

  expect_s3_class(bins, "data.frame")
  expect_true(nrow(bins) > 0)
  expect_false(any(is.na(bins$chr)))
})

test_that("CountRawReads rejects invalid bin sizes before reading BAM headers", {
  bam_file <- tempfile(fileext = ".bam")
  file.create(bam_file)

  expect_error(CountRawReads(bam_file, binSize = 199), "Recommended binSize")
  expect_error(CountRawReads(bam_file, binSize = 10001), "Recommended binSize")
})

test_that("ReadMeta rejects any duplicated sample ID", {
  meta_file <- tempfile(fileext = ".txt")
  writeLines(
    c(
      "ID\tANTIBODY\tGROUP",
      "sample.bam\tH3K27me3\tWT",
      "sample.bam\tH3K27me3\tK27M"
    ),
    meta_file
  )

  expect_error(ReadMeta(meta_file), "duplicate ID")
})

test_that("CalculateSF chooses scaling reference from passing samples only", {
  cutoff <- seq(0, 10, by = 0.1)
  parsed <- data.frame(
    cutoff = cutoff,
    pass_fast = pmin(cutoff / 4, 0.995),
    pass_slow = pmin(cutoff / 8, 0.995),
    fail_high_slope = pmin(cutoff / 0.7, 0.995),
    check.names = FALSE
  )
  meta <- data.frame(
    ID = c("pass_fast", "pass_slow", "fail_high_slope"),
    ANTIBODY = "H3K27me3",
    GROUP = c("WT", "K27M", "INPUT"),
    COLOR = c("grey", "green", "black"),
    stringsAsFactors = FALSE
  )
  rownames(meta) <- meta$ID

  sf <- CalculateSF(
    data = parsed,
    metaFile = meta,
    minFirstTurn = 0.5,
    maxLastTurn = 0.99,
    cutoff_QC = 1.2
  )
  passing <- sf$QC == "pass"

  expect_true(any(passing))
  expect_true(any(!passing))
  expect_equal(min(sf$SF[passing], na.rm = TRUE), 1)
  expect_true(all(is.na(sf$SF[!passing])))
})

test_that("ChIPseqSpikeInFree invisibly returns calculated scaling factors", {
  meta <- data.frame(
    ID = "sample.bam",
    ANTIBODY = "H3K27me3",
    GROUP = "WT",
    COLOR = "grey",
    stringsAsFactors = FALSE
  )
  rownames(meta) <- meta$ID
  raw_counts <- data.frame(sample.bam = c(1000001, 1000002), check.names = FALSE)
  parsed <- data.frame(cutoff = c(0, 1), sample.bam = c(0, 0.95), check.names = FALSE)
  sf <- cbind(meta, QC = "pass", TURNS = "0.5,0.2,1.5,0.95", SF = 1)

  fn_env <- environment(ChIPseqSpikeInFree)
  mocked <- list(
    ReadMeta = function(metaFile) meta,
    CountRawReads = function(...) raw_counts,
    ParseReadCounts = function(...) parsed,
    CalculateSF = function(...) sf,
    PlotDistr = function(...) NULL,
    BoxplotSF = function(...) NULL
  )
  old <- mget(names(mocked), envir = fn_env, inherits = FALSE)
  on.exit({
    for (nm in names(old)) {
      assign(nm, old[[nm]], envir = fn_env)
    }
  }, add = TRUE)
  for (nm in names(mocked)) {
    assign(nm, mocked[[nm]], envir = fn_env)
  }

  result <- ChIPseqSpikeInFree(
    bamFiles = "sample.bam",
    chromFile = "hg19",
    metaFile = "sample_meta.txt",
    prefix = tempfile(),
    disableOverwritten = FALSE
  )

  expect_equal(result, sf)
})

test_that("ChIPseqSpikeInFree rejects invalid bin size early", {
  expect_error(
    ChIPseqSpikeInFree(bamFiles = character(), binSize = 199),
    "recommended binSize"
  )
})