if (!exists("GenerateBins", mode = "function")) {
  source_file <- "R/ChIPseqSpikeInFree.R"
  if (!file.exists(source_file)) {
    source_file <- "../../R/ChIPseqSpikeInFree.R"
  }
  source(source_file, local = globalenv())
}