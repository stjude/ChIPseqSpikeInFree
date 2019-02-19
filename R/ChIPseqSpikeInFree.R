######################################################
#' generate genome-wide bins for counting purpose
#'
#' Given a chrom.size file, this function allows you to generate a your own sliding windows (bins).
#' @param chromFile chrom.size. Given "hg19","mm10","mm9" or "hg38", will load chrom.size file from package folder.
#' @param binSize size of bins (bp)
#' @param overlap overlaps between two consecutive bins (bp)
#' @param withChr chromosome names in bin File have chr if set withChr to TRUE; FALSE - no chr
#' @return A data.frame of generated bins
#' @export
#' @examples
#' ## 1. generate a mm10 binFile without chr and use a binSize of 1000 bp
#' ## and overlap of 500 bp between two consecutive bins
#' 
#' ## "mm10" will be parsed as system.file("extdata", "mm10.chrom.sizes", 
#' ##       package = "ChIPseqSpikeInFree")
#' # binDF <- GenerateBins(chromFile="mm10",binSize=1000, overlap=500,
#' # withChr=FALSE,prefix="mm10")
#' 
#' ## 2. generate a hg19 binFile with chr and use a binSize of 2000 bp
#' 
#' ## "hg19" will be parsed as system.file("extdata", "hg19.chrom.sizes",
#' ##       package = "ChIPseqSpikeInFree")
#' # binDF<- GenerateBins(chromFile="hg19",binSize=2000, overlap=0,
#' # withChr=TRUE,prefix="hg19")
GenerateBins <- function(chromFile, binSize = 1000, overlap = 0, withChr = TRUE) {
  # generate genome-wide bins for counting purpose

  options(stringsAsFactors = FALSE)
  if (!file.exists(chromFile)) {
    chromFile <- system.file("extdata", paste0(chromFile, ".chrom.sizes"), package = "ChIPseqSpikeInFree")
    if (!file.exists(chromFile)) {
      stop(paste0("\nchromFile was not found:", chromFile, "\n"))
    }
  }
  cat(paste0("\n\tFollowing file will be used to generate sliding windows: ", basename(chromFile)))
  chromSize <- read.table(chromFile, sep = "\t", header = F, fill = TRUE, quote = "")
  options(scipen = 999) # disable scientific notation when writing out
  counter <- NULL
  excluded <- grepl("random", chromSize[, 1]) | grepl("Un", chromSize[, 1]) | grepl("alt", chromSize[, 1]) | grepl("hap", chromSize[, 1])
  chromSize <- chromSize[!excluded, ]
  chrNotation <- gsub("chrM", "MT", chromSize[, 1])
  chrNotation <- gsub("chr", "", chrNotation)
  if (withChr) {
    chrNotation <- gsub("^", "chr", chrNotation)
    chrNotation <- gsub("chrMT", "chrM", chrNotation)
  }
  chromSize[, 1] <- chrNotation
  beds <- NULL
  for (r in 1:nrow(chromSize)) {
    chrName <- chromSize[r, 1]
    chrLength <- chromSize[r, 2]
    realBinSize <- binSize - overlap
    starts <- seq(0, chrLength, realBinSize) + 1
    ends <- c(seq(binSize, chrLength, realBinSize), chrLength)
    starts <- starts[1:length(ends)]
    bed <- data.frame(chr = chrName, start = starts, end = ends)
    counter <- c(counter, nrow(bed))
    if (r == 1) { # overwrite lines
      beds <- bed
    } else {
      beds <- rbind(beds, bed)
    }
  }
  return(beds)
}

######################################################
#' count raw reads for all bins
#'
#' This function counts raw reads for each bin.
#' @param bamFiles a vector of bam filenames.
#' @param chromFile a chrom.size file. Given "hg19","mm10","mm9" or "hg38", will load chrom.size file from package folder. Otherwise, give a /your/path/chrom.size
#' @param prefix prefix of output file name
#' @param singleEnd To count paired-end reads, set argument singleEnd=FALSE
#' @param binSize size of bins (bp)
#' @return a data.frame of raw counts for each bin
#' @import Rsamtools
#' @import GenomicRanges
#' @import GenomicAlignments
#' @import IRanges
#' @export
#' @examples
#' ## 1.count reads using mm9 bams
#' # bams <- c("your/path/ChIPseq1.bam","your/path/ChIPseq2.bam")
#' # rawCountDF <- CountRawReads(bamFiles=bams,chromFile="mm9",prefix="your/path/test",singleEnd=TRUE)
CountRawReads <- function(bamFiles, chromFile = "hg19", prefix = "test", singleEnd = TRUE, binSize=1000) {
  # count raw reads for every 1kb bin across genome

  # check chromosome notation in bam file
  bamHeader <- scanBamHeader(bamFiles[1])
  chrFlag <- grepl("SN:chr", bamHeader[[1]]$text[2])
  binDF <- GenerateBins(chromFile = chromFile, binSize = 1000, overlap = 0, withChr = chrFlag)
  bins <- binDF
  myCoords <- GRanges(
    seqnames = bins[, 1],
    ranges = IRanges(start = bins[, 2], end = bins[, 3], names = paste0(bins[, 1], ":", bins[, 2], "-", bins[, 3])),
    strand = "+"
  )

  cat("\n\tThis step could take some time. Please be patient...")
  bamlist <- BamFileList(bamFiles, yieldSize = 5e4)
  names(bamlist) <- basename(bamFiles)
  assay <- NULL ## To please R CMD check
  counts <- summarizeOverlaps(
    features = myCoords,
    reads = bamlist,
    ignore.strand = TRUE,
    singleEnd = singleEnd
  )
  dat <- as.data.frame(assay(counts))
  failedBams <- sum(colSums(dat) == 0)
  if (failedBams > 0) {
    cat("\n\tWarning:  read counts are all zeros in ", failedBams, " column(s). \n")
  }
  datOut <- data.frame(bin = rownames(dat), dat, check.names = F)
  options(scipen = 99) # disable scientific notation when writing out
  outFile <- paste0(prefix, "_rawCounts.txt")
  write.table(datOut, outFile, sep = "\t", quote = F, row.names = F, col.names = T)
  invisible(return(dat))
}


######################################################
#' read in sample metadata file
#'
#' This function allows you to load metadat to a R data.frame and return the object.
#' In addtion, it validates meta_info format and adds a COLOR column if it's undefined.
#' @param metaFile  a metadata file name; the file must have three columns: ID (bam filename without full path), ANTIBODY and GROUP. the COLOR column is optional and will be used for plotting purpose.
#' @return A data.frame of metaFile
#' @export
#' @examples
#' ## 1. load an example of metadata file
#' 
#' metaFile <- system.file("extdata", "sample_meta.txt", package = "ChIPseqSpikeInFree")
#' meta <- ReadMeta(metaFile)
#' head(meta, n = 1)
#' meta
#' #                                             ID ANTIBODY GROUP COLOR
#' #  H3K27me3-NSH.K27M.A.bam H3K27me3-NSH.K27M.A.bam H3K27me3  K27M   red
#' #  H3K27me3-NSH.K27M.B.bam H3K27me3-NSH.K27M.B.bam H3K27me3  K27M   red
#' #  H3K27me3-NSH.K27M.C.bam H3K27me3-NSH.K27M.C.bam H3K27me3  K27M   red
#' #  H3K27me3-NSH.WT.D.bam     H3K27me3-NSH.WT.D.bam H3K27me3    WT  blue
#' #  H3K27me3-NSH.WT.E.bam     H3K27me3-NSH.WT.E.bam H3K27me3    WT  blue
#' #  H3K27me3-NSH.WT.F.bam     H3K27me3-NSH.WT.F.bam H3K27me3    WT  blue
ReadMeta <- function(metaFile = "sample_meta.txt") {
  # read in sample metadata file

  if (!file.exists(metaFile)) {
    stop(paste0("\nFile was not found:", metaFile, "\n"))
  }
  meta <- read.table(metaFile, sep = "\t", header = TRUE, fill = TRUE, stringsAsFactors = FALSE, quote = "", check.names = F)
  colnames(meta) <- toupper(colnames(meta))
  rownames(meta) <- meta$ID
  if (!"ID" %in% colnames(meta) | !"ANTIBODY" %in% colnames(meta) | !"GROUP" %in% colnames(meta)) {
    cat("\n**ERROR: Invalid meta file. **\n\tID, ANTIBODY, GROUP column(s) are required. COLOR column is optional.\n")
    # stop("Exit due to unexpected error.")
    return(NA)
  }
  if (!"COLOR" %in% colnames(meta)) {
    myColors <- c(
      "red", "blue", "green", "#F4B5BE", "#5F5647",
      "#9C974A", "#9B110E", "#CDBF8B", "#3F5151",
      "#F9D77A", "#4E2A1E", "#85D4E3", "#0C1707", "yellow", "grey"
    )
    groupNum <- as.numeric(factor(paste(meta$ANTIBODY, meta$GROUP, sep = "")))
    meta$COLOR <- myColors[groupNum]
  }
  invisible(meta)
}

######################################################

#' parse readCounts matrix
#'
#' This function allows you to parse rawCount table (generated by CountRawReads() function) to a parsedMatrix of (cutoff, and percent of reads accumulatively passed the cutoff in each sample).
#' @param data a data.frame returned by readRawCounts() or a file name of rawCount table
#' @param metaFile a data.frame of metadata by ReadMeta(); or a filename of metadata file.
#' @param by step used to define cutoffs; ParseReadCounts will cumulatively calculate the percent of reads that pass the every cutoff.
#' @param prefix prefix of output filename to save the parsedMatrix of (cutoff, and percent of reads accumulatively passed the cutoff in each sample).
#' @return A data.frame of parsed data.
#' @export
#' @examples
#' 
#' ## prerequisite step 1. count raw reads
#' ## (if your bam files were aligned to mm9 genome with chr in reference chromosomes).
#' # bams <- c("your/path/ChIPseq1.bam","your/path/ChIPseq2.bam")
#' # rawCountDF <- CountRawReads(bamFiles=bams,chromFile="mm9",prefix="your/path/test")
#' ## output file will be "your/path/test_rawCount.txt"
#' # head(rawCountDF,n=2)
#' 
#' # bin   ChIPseq1.bam    ChIPseq2.bam
#' # chr1:1-1000  0   0
#' # chr1:1001-2000   0   0
#' 
#' ## prerequisite step 2: generate your sample_meta.txt.
#' ##  A tab-delimited txt file has three required columns
#' # ID    ANTIBODY    GROUP
#' # ChIPseq1.bam  H3K27me3    WT
#' # ChIPseq2.bam  H3K27me3    K27M
#' 
#' ## 1.parse readCount table using this function.
#' # metaFile <- "your/path/sample_meta.txt"
#' # dat <- ParseReadCounts(data="your/path/test_rawCount.txt",
#' # metaFile=metaFile, prefix="your/path/test")
#' ## output file will be "your/path/test_parsedMatrix.txt"
ParseReadCounts <- function(data, metaFile = "sample_meta.txt", by = 0.05, prefix = "test") {
  if (class(metaFile) == "character") { # given a filename, need to load it
    metaFile <- ReadMeta(metaFile)
  }
  if (class(data) != "data.frame") {
    if (file.exists(data)) {
      data <- read.table(data, sep = "\t", header = TRUE, fill = TRUE, stringsAsFactors = FALSE, quote = "", row.names = 1, check.names = F)
    } else {
      stop("Not found rawCount table and your input is not a data.frame.\n")
    }
  }
  kept <- intersect(metaFile$ID, colnames(data))
  if (length(kept) == 0) {
    cat("\n** Oooops: you need to change metadata file **\n")
    stop("Please check whether you IDs in metaFile match with colnames(bam filenames) in parsedMatrix.\n\n")
  }
  metaFile <- metaFile[kept, ] # only use samples listed in metaFile
  CPM <- apply(data, 2, function(x) {
    x / sum(x) * 1000000
  })
  colnames(CPM) <- colnames(data)
  MAX <- apply(CPM, 2, max)
  MAX <- ceiling(max(MAX))
  SEQ <- seq(0, MAX, by = by) # smaller value, higher resolution
  dat <- data.frame(cutoff = SEQ)
  for (ind in 1:ncol(CPM)) {
    cat("\n", colnames(CPM)[ind], "\n\t")
    tmp <- CPM[, ind]
    totalCPM <- sum(tmp)
    res <- NULL
    counter <- 0
    for (K in SEQ) {
      counter <- counter + 1
      if (counter %% 100 == 0 | counter == 1) {
        cat(paste0(counter, ".."))
      }
      pct <- sum(tmp[tmp <= K]) / totalCPM
      res <- c(res, pct)
      if (pct == 1) { # exit loop when all reads have been used
        res <- c(res, rep(1, nrow(dat) - length(res)))
        break
      }
    }
    dat <- cbind(dat, res)
    colnames(dat)[ncol(dat)] <- colnames(CPM)[ind]
  }
  output <- paste0(prefix, "_parsedMatrix.txt")
  write.table(dat, output, sep = "\t", quote = F, row.names = F, col.names = T)
  cat("\n\t", output, "[saved]")

  return(dat)
}


#' calculate scaling factors and save results
#'
#' This function allows you to plot curves, caculate SF(scaling factor) per antibody based on parsedMatrix.
#' In addtion return a data.frame of updated metadata.
#' If you run ChIPseqSpikeInFree() seperately for two batches, the scaling factors will be not comparable between two batches.
#' The correct way is to combine bamFiles parameter and create a new metadata file to include all bam files. Then re-run ChIPseqSpikeInFree().
#' @param data a data.frame generated by function ParseReadCounts() or a file name of parsed matrix
#' @param metaFile a data.frame of metadata by ReadMeta(); or a filename of metadata file.
#' @param dataRange set c(0.1 and 0.99) by default. The first 10 percent of noisy weakest reads were excluded and stop calculation when 99 percent of reads are used.
#' @param prefix prefix of output filename to save the scaling factor values.
#' @return A data.frame of the updated metaFile with scaling factors
#' @export
#' @examples
#' ## 1. start from a parsedMatrix file
#' 
#' # parsedMatrixFile <- "your/path/test_parsedMatrix.txt"
#' # metaFile <- "your/path/sample_meta.txt"
#' # parsedDF <- read.table(parsedMatrixFile, sep="\t",header=TRUE,fill=TRUE,
#' # quote="",row.names=NULL ,check.names=F)
#' # res <- CalculateSF (data=parsedDF,metaFile=metaFile, prefix="your/path/test")
#' 
#' ## 2. start from a rawCount file
#' 
#' # metaFile <- "your/path/sample_meta.txt"
#' # parsedDF <- ParseReadCounts(data="your/path/test_rawCounts.txt", metaFile=metaFile,
#' #     prefix="your/path/test_parsedMatrix.txt")
#' # res <- CalculateSF (data=parsedDF,metaFile=metaFile,
#' #     prefix="your/path/test")
CalculateSF <- function(data, metaFile = "sample_meta.txt", prefix = "test", dataRange = c(0.1, 0.99)) {
  # calculate scaling factors

  pctMin <- dataRange[1]
  pctMAX <- dataRange[2]
  if (class(metaFile) == "character") { # given a filename, need to load it
    meta <- ReadMeta(metaFile)
  } else {
    meta <- metaFile
  }
  if (class(data) == "character") { # given a filename, need to load it
    if (file.exists(data)) {
      data <- read.table(data, sep = "\t", header = TRUE, fill = TRUE, stringsAsFactors = FALSE, quote = "", row.names = NULL, check.names = F)
    } else {
      stop("Not found parsed matrix and your input is not a data.frame.\n")
    }
  }
  kept <- intersect(meta$ID, colnames(data))
  if (length(kept) == 0) {
    stop("meta data IDs don't match with column names in parsedMatrix.txt.")
  }
  meta <- meta[kept, ] # only use samples listed in meta
  data <- data[, c(colnames(data)[1], kept)] # reorder or subset
  imgOutput <- paste0(prefix, "_distribution.pdf")
  pdf(imgOutput, width = 6, height = 6)
  slopes <- NULL
  MAX_CPM <- max(data[,1])
  for (r in 2:ncol(data)) {
    x <- data[, 1]
    y <- data[, r]
    used <- data.frame(x = x, y = y)
    used <- na.omit(used)
    yMin <- y[sum(y < pctMin)]
    yMax <- y[sum(y < pctMAX)]
    xMin <- min(used$x[used$y == yMin])
    xMax <- max(used$x[used$y == yMax])
    slope <- round((yMax - yMin) / (xMax - xMin), 5)
    slopes <- c(slopes, slope)
    if (r == 2) {
      plot(used,
        main = imgOutput, col = meta$COLOR[r - 1], lwd = 2, xlab = "cutoff (CPMW)",
        cex = 0.8, cex.main = 0.5, type = "l", ylab = "Proportion of reads", xlim = c(0, MAX_CPM), ylim = c(0.0, 1.1)
      )
    } else {
      lines(used, col = meta$COLOR[r - 1], lty = 1, lwd = 2, pch = 20, cex = 0.1)
    }
    lines(x = c(xMin, xMax), y = c(yMin, yMax), col = "blue", lty = 3)
  }
  meta$SLOPE <- slopes
  meta$SF <- NA
  for (ab in unique(meta$ANTIBODY)) {
    inds <- grep(ab, meta$ANTIBODY)
    SF <- round(max(slopes[inds]) / slopes[inds], 2)
    meta$SF[inds] <- SF
  }

  legend("bottomright", legend = paste(meta$ID, paste(",sf=", meta$SF, sep = ""), sep = ""), col = meta$COLOR, pch = 15, bty = "n", ncol = 1, cex = 0.8)
  garbage <- dev.off()
  cat("\n\t", imgOutput, "[saved]")
  output <- paste0(prefix, "_SF.txt")
  outDF <- meta[, !colnames(meta) %in% "SLOPE"]
  write.table(outDF, output, sep = "\t", quote = F, row.names = F, col.names = T)
  cat("\n\t", output, "[saved]")
  cat("\n\n  Reporting summary")
  if ("GROUP" %in% colnames(meta)) {
    combs <- paste(meta$ANTIBODY, ".", meta$GROUP, sep = "")
    groups <- unique(combs)
    for (group in groups) {
      inds <- grep(group, combs)
      cat("\n\t", group, ",\tave.SF =", round(mean(meta$SF[inds]), 3))
      #  cat(",\tave.SLOPE =", round(mean(meta$SLOPE[inds]),4))
    }
    cat("\n")
  }

  invisible(return(meta))
}
######################################################

#' This function generates boxplot using sacaling factor table. It's been included in the last step of ChIPseqSpikeInFree().
#'
#' @param input a file/data.frame generated/returned by CalculateSF() or ChIPseqSpikeInFree().
#' It looks like metadata file but has extra columns SF and COLOR.
#' @param prefix prefix of output filename.
#' @return A filename of generated boxplot
#' @export
#' @examples
#' ## 1. re-generate boxplot of ChIPseqSpikeInFree scaling factors
#' ## After you run ChIPseqSipkeFree(), a sacaling factor table
#' ## (for example, test_SF.txt) will be generated.
#' 
#' # BoxplotSF(input="test_SF.txt",prefix="test")
BoxplotSF <- function(input, prefix = "test") {
  # This function generates boxplot using sacaling factor table
  if (class(input) == "character") {
    input <- read.table(input, sep = "\t", header = TRUE, fill = TRUE, stringsAsFactors = FALSE, quote = "", row.names = 1, check.names = F)
  }

  if (!"SF" %in% colnames(input) || !"ANTIBODY" %in% colnames(input) || !"GROUP" %in% colnames(input)) {
    stop("Input looks invalid for BoxplotSF().\n")
  }
    
  input$GROUP2 <- paste(input$GROUP,input$ANTIBODY,sep=".")
  input <- input[order(input$ANTIBODY,input$GROUP), ]
  groupLabels <- unique(input$GROUP2)
  if (!"COLOR" %in% colnames(input)) {
    tim10equal <- c("skyblue", "#EF0000", "grey", "#00DFFF", "#50FFAF", "#BFFF40", "#FFCF00", "#FF6000", "#0000FF", "#800000")
    myCols <- tim10equal[as.factor(input$GROUP2)]
  } else {
    myCols <- input[!duplicated(input$GROUP2), "COLOR"]
  }
  output <- paste0(prefix, "_boxplot.pdf")
  pdf(output, width = 7, height = 6)
  par(mfrow = c(1, 1), oma = c(0, 0, 0, 0), mar = c(10, 4, 2, 12))
  myTitle <- "ScalingFactor~GROUP+ANTIBODY"
  nameOrder <- ordered(input$GROUP2, levels = groupLabels)
  bp <- boxplot(SF ~ nameOrder,
    data = input, main = myTitle,
    las = 2, ylab = "Scaling Factors", xlab = "", cex.axis = 0.8, outline = T, col = myCols,
    cex.lab = 0.7, cex.main = 0.8
  )
  stripchart(SF ~ nameOrder,
    data = input[!input$SF %in% bp$out, ], vertical = T, 
    method = "jitter", add = TRUE, jitter = 0.2, cex = 0.7, pch = 1, col = "#595959"
  )
  par(xpd = T)
  legend(x = length(groupLabels) + 1, y = max(input$SF), legend = groupLabels, col = myCols, pch = 15, bty = "n", ncol = 1, cex = )
  garbage <- dev.off()
  cat("\n\t", output, "[saved]")
  invisible(output)
}

######################################################

#' wrapper function - perform ChIP-seq spike-free normalization in one step.
#'
#' This function wraps all steps.
#' If you run ChIPseqSpikeInFree() seperately for two batches, the scaling factors will be not comparable between two batches.
#' The correct way is to combine bamFiles parameter and create a new metadata file to include all bam files. Then re-run ChIPseqSpikeInFree().
#' @param bamFiles a vector of bam filenames.
#' @param chromFile chrom.size file. Given "hg19","mm10","mm9" or "hg38", will load chrom.size file from package folder.
#' @param metaFile a filename of metadata file. the file must have three columns: ID (bam filename without full path), ANTIBODY and GROUP
#' @param binSize size of bins (bp)
#' @param prefix prefix of output filename.
#' @return A data.frame of the updated metaFile with scaling factor
#' @export
#' @examples
#' ## 1 first You need to generate a sample_meta.txt (tab-delimited txt file).
#' # metaFile <- "your/path/sample_meta.txt"
#' # meta <- ReadMeta(metaFile)
#' # head(meta)
#' # ID ANTIBODY GROUP
#' # ChIPseq1.bam H3K27me3 WT
#' # ChIPseq2.bam H3K27me3 K27M
#' 
#' ## 2. bam files
#' # bams <- c("ChIPseq1.bam","ChIPseq2.bam")
#' # prefix <- "test"
#' 
#' ## 3. run ChIPseqSpikeInFree pipeline
#' # ChIPseqSpikeInFree(bamFiles=bams, chromFile="mm9",metaFile=metaFile,prefix="test")
ChIPseqSpikeInFree <- function(bamFiles, chromFile = "hg19", metaFile = "sample_meta.txt", prefix = "test", binSize=1000) {
  # perform ChIP-seq spike-free normalization in one step

  cat("\nstep1. loading metadata file...")
  meta <- ReadMeta(metaFile)
  cat("\n\t[--done--]\n")
  cat("\nstep2. counting reads...")
  output1 <- paste0(prefix, "_rawCounts.txt")
  disableOverwritten <- TRUE

  # in case of downstream errors, just load existing outFile to avoid re-counting
  if (file.exists(output1) && disableOverwritten) {
    cat("\n\t", output1, "[just loading the existing file; delete this file first to do re-counting ]")
    rawCountDF <- read.table(output1, sep = "\t", header = TRUE, fill = TRUE, stringsAsFactors = FALSE, quote = "", row.names = 1, check.names = F)
  } else {
    rawCountDF <- CountRawReads(bamFiles = bamFiles, chromFile = chromFile, prefix = prefix, binSize=binSize)
  }
  cat("\n\t[--done--]\n")
  cat("\nstep3. parsing raw counts...")
  output2 <- paste0(prefix, "_parsedMatrix.txt")

  # in case of downstream errors, just load existing outFile to avoid re-processing raw counts
  if (file.exists(output2) && disableOverwritten) {
    cat("\n\t", output1, "[just loading the existing file; delete this file first to re-parse rawCount table]")
    parsedDF <- read.table(output2, sep = "\t", header = TRUE, fill = TRUE, stringsAsFactors = FALSE, quote = "", row.names = NULL, check.names = F)
  } else {
    steps <- 0.05* binSize/1000
    parsedDF <- ParseReadCounts(data = rawCountDF, metaFile = meta, prefix = prefix, by = steps)
  }
  cat("\n\t[--done--]\n")
  cat("\nstep4. calculating scaling factors...")
  result <- CalculateSF(data = parsedDF, metaFile = meta, prefix = prefix, dataRange = c(0.1, 0.99))
  cat("\t[--done--]\n")
  cat("\nstep5. ploting scaling factors...")
  BoxplotSF(result, prefix)
  cat("\n\t[--done--]\n")
  invisible(result)
}
