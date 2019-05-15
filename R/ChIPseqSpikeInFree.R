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
  if (binSize < 200 && binSize > 10000) {
    stop(paste0("\n**Recommended binSize range 200 ~ 10000 (bp); Your binSize is", binSize, " **\n"))
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
#' @param binSize size of bins (bp). Recommend a value bwteen 200 and 10000
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
CountRawReads <- function(bamFiles, chromFile = "hg19", prefix = "test", singleEnd = TRUE, binSize = 1000) {
  # count raw reads for every 1kb bin across genome
  # check file
  if (sum(!file.exists(bamFiles)) > 0) {
    cat("\n** some bam files are unavailable:**")
    cat("\t", bamFiles[!file.exists(bamFiles)], "\n")
    stop("Invalid bam file list!")
  }

  if (binSize < 200 && binSize > 10000) {
    stop(paste0("\n**Recommended binSize range 200 ~ 10000 (bp); Your binsize is", binSize, " **\n"))
  }
  # Function to check chromosome notation in bam file
  # return a logical vector
  ChrCheck <- function(bamFiles) {
    bamHeader <- scanBamHeader(bamFiles)
    nameList <- unlist(
      lapply(
        1:length(bamHeader),
        function(x) {
          names(bamHeader[[x]]$targets[1])
        }
      )
    )
    grepl("chr", nameList)
  }

  # Function to count reads
  # return a data frame
  CountingReads <- function(bamFiles, withChr) {
    binDF <- GenerateBins(chromFile = chromFile, binSize = binSize, overlap = 0, withChr = withChr)
    bins <- binDF
    myCoords <- GRanges(
      seqnames = bins[, 1],
      ranges = IRanges(start = bins[, 2], end = bins[, 3], names = paste0(bins[, 1], ":", bins[, 2], "-", bins[, 3])),
      strand = "+"
    )
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
    return(dat)
  }

  chrFlags <- ChrCheck(bamFiles)
  bams.Chr <- bamFiles[chrFlags]
  bams.NoChr <- bamFiles[!chrFlags]

  cat("\n\tThis step could take some time. Please be patient...")
  dat.Chr <- NULL
  dat.NoChr <- NULL
  if (length(bams.Chr) > 0) {
    dat.Chr <- CountingReads(bams.Chr, TRUE)
  }
  if (length(bams.NoChr) > 0) {
    dat.NoChr <- CountingReads(bams.NoChr, FALSE)
  }
  if (!is.null(dat.Chr) && !is.null(dat.NoChr)) {
    datOut <- data.frame(bin = rownames(dat.Chr), dat.Chr, dat.NoChr, check.names = F)
  } else if (!is.null(dat.Chr)) {
    datOut <- data.frame(bin = rownames(dat.Chr), dat.Chr, check.names = F)
  } else {
    datOut <- data.frame(bin = rownames(dat.NoChr), dat.NoChr, check.names = F)
  }
  rm(list = c("dat.Chr", "dat.NoChr"))
  options(scipen = 99) # disable scientific notation when writing out
  outFile <- paste0(prefix, "_rawCounts.txt")
  write.table(datOut, outFile, sep = "\t", quote = F, row.names = F, col.names = T)
  invisible(return(datOut[, -1]))
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
  # check duplicate ID
  if (sum(duplicated(meta$ID)) > 1) {
    cat("\n**Warning: found duplicate ID(s) in your metaFile.**\n\n")
    print(unique(meta$ID[duplicated(meta$ID)]))
    cat("\nPlease correct your metaFile and then re-run the pipeline.\n\n")
    quit("no")
  }
  rownames(meta) <- meta$ID
  if (!"ID" %in% colnames(meta) | !"ANTIBODY" %in% colnames(meta) | !"GROUP" %in% colnames(meta)) {
    cat("\n**ERROR: Invalid meta file. **\n\tID, ANTIBODY, GROUP column(s) are required. COLOR column is optional.\n")
    # stop("Exit due to unexpected error.")
    return(NA)
  }
  if (!"COLOR" %in% colnames(meta)) {
    # myColors <- c("red", "blue", "green", "#F4B5BE", "#5F5647", "#9C974A", "#9B110E", "#CDBF8B", "#3F5151","#F9D77A", "#4E2A1E", "#85D4E3", "#0C1707", "yellow", "grey" )
    myColors <- c(
      "turquoise", "blue", "brown", "yellow", "green", "red", "black", "pink", "magenta", "purple", "greenyellow", "tan", "salmon", "cyan", "midnightblue", "lightcyan",
      "grey60", "lightgreen", "lightyellow", "royalblue", "darkred", "darkgreen", "darkturquoise", "darkgrey", "orange", "darkorange", "white", "skyblue",
      "saddlebrown", "steelblue", "paleturquoise", "violet", "darkolivegreen", "darkmagenta", "sienna3", "yellowgreen", "skyblue3", "plum1", "orangered4", "mediumpurple3",
      "lightsteelblue1", "lightcyan1", "ivory", "floralwhite", "darkorange2", "brown4", "bisque4", "darkslateblue", "plum2", "thistle2", "thistle1", "salmon4",
      "palevioletred3", "navajowhite2", "maroon", "lightpink4", "lavenderblush3", "honeydew1", "darkseagreen4", "coral1", "antiquewhite4", "coral2", "mediumorchid", "skyblue2",
      "yellow4", "skyblue1", "plum", "orangered3", "mediumpurple2", "lightsteelblue", "lightcoral", "indianred4", "firebrick4", "darkolivegreen4", "brown2", "blue2",
      "darkviolet", "plum3", "thistle3", "thistle", "salmon2", "palevioletred2", "navajowhite1", "magenta4", "lightpink3", "lavenderblush2", "honeydew", "darkseagreen3",
      "coral", "antiquewhite2", "coral3", "mediumpurple4", "skyblue4", "yellow3", "sienna4", "pink4", "orangered1", "mediumpurple1", "lightslateblue", "lightblue4",
      "indianred3", "firebrick3", "darkolivegreen2", "blueviolet", "blue4", "deeppink", "plum4", "thistle4", "tan4", "salmon1", "palevioletred1", "navajowhite",
      "magenta3", "lightpink2", "lavenderblush1", "green4", "darkseagreen2", "chocolate4", "antiquewhite1", "coral4", "mistyrose", "slateblue", "yellow2", "sienna2",
      "pink3", "orangered", "mediumpurple", "lightskyblue4", "lightblue3", "indianred2", "firebrick2", "darkolivegreen1", "blue3", "brown1", "deeppink1", "powderblue",
      "tomato", "tan3", "royalblue3", "palevioletred", "moccasin", "magenta2", "lightpink1", "lavenderblush", "green3", "darkseagreen1", "chocolate3", "aliceblue",
      "cornflowerblue", "navajowhite3", "slateblue1", "whitesmoke", "sienna1", "pink2", "orange4", "mediumorchid4", "lightskyblue3", "lightblue2", "indianred1", "firebrick",
      "darkgoldenrod4", "blue1", "brown3", "deeppink2", "purple2", "tomato2", "tan2", "royalblue2", "paleturquoise4", "mistyrose4", "magenta1", "lightpink",
      "lavender", "green2", "darkseagreen", "chocolate2", "antiquewhite", "cornsilk", "navajowhite4", "slateblue2", "wheat3", "sienna", "pink1", "orange3",
      "mediumorchid3", "lightskyblue2", "lightblue1", "indianred", "dodgerblue4", "darkgoldenrod3", "blanchedalmond", "burlywood", "deepskyblue", "red1", "tomato4", "tan1",
      "rosybrown4", "paleturquoise3", "mistyrose3", "linen", "lightgoldenrodyellow", "khaki4", "green1", "darksalmon", "chocolate1", "antiquewhite3", "cornsilk2", "oldlace",
      "slateblue3", "wheat1", "seashell4", "peru", "orange2", "mediumorchid2", "lightskyblue1", "lightblue", "hotpink4", "dodgerblue3", "darkgoldenrod1", "bisque3",
      "burlywood1", "deepskyblue4", "red4", "turquoise2", "steelblue4", "rosybrown3", "paleturquoise1", "mistyrose2", "limegreen", "lightgoldenrod4", "khaki3", "goldenrod4",
      "darkorchid4", "chocolate", "aquamarine", "cyan1", "orange1", "slateblue4", "violetred4", "seashell3", "peachpuff4", "olivedrab4", "mediumorchid1", "lightskyblue",
      "lemonchiffon4", "hotpink3", "dodgerblue1", "darkgoldenrod", "bisque2", "burlywood2", "dodgerblue2", "rosybrown2", "turquoise4", "steelblue3", "rosybrown1", "palegreen4",
      "mistyrose1", "lightyellow4", "lightgoldenrod3", "khaki2", "goldenrod3", "darkorchid3", "chartreuse4", "aquamarine1", "cyan4", "orangered2", "snow", "violetred2",
      "seashell2", "peachpuff3", "olivedrab3", "mediumblue", "lightseagreen", "lemonchiffon3", "hotpink2", "dodgerblue", "darkblue", "bisque1", "burlywood3", "firebrick1",
      "royalblue1", "violetred1", "steelblue1", "rosybrown", "palegreen3", "mintcream", "lightyellow3", "lightgoldenrod2", "khaki1", "goldenrod2", "darkorchid2", "chartreuse3",
      "aquamarine2", "darkcyan", "orchid", "snow2", "violetred", "seashell1", "peachpuff2", "olivedrab2", "mediumaquamarine", "lightsalmon4", "lemonchiffon2", "hotpink1",
      "deepskyblue3", "cyan3", "bisque", "burlywood4", "forestgreen", "royalblue4", "violetred3", "springgreen3", "red3", "palegreen1", "mediumvioletred", "lightyellow2",
      "lightgoldenrod1", "khaki", "goldenrod1", "darkorchid1", "chartreuse2", "aquamarine3", "darkgoldenrod2", "orchid1", "snow4", "turquoise3", "seashell", "peachpuff1",
      "olivedrab1", "maroon4", "lightsalmon3", "lemonchiffon1", "hotpink", "deepskyblue2", "cyan2", "beige", "cadetblue", "gainsboro", "salmon3", "wheat",
      "springgreen2", "red2", "palegreen", "mediumturquoise", "lightyellow1", "lightgoldenrod", "ivory4", "goldenrod", "darkorchid", "chartreuse1", "aquamarine4", "darkkhaki",
      "orchid3", "springgreen1", "turquoise1", "seagreen4", "peachpuff", "olivedrab", "maroon3", "lightsalmon2", "lemonchiffon", "honeydew4", "deepskyblue1", "cornsilk4",
      "azure4", "cadetblue1", "ghostwhite", "sandybrown", "wheat2", "springgreen", "purple4", "palegoldenrod", "mediumspringgreen", "lightsteelblue4", "lightcyan4", "ivory3",
      "gold3", "darkorange4", "chartreuse", "azure", "darkolivegreen3", "palegreen2", "springgreen4", "tomato3", "seagreen3", "papayawhip", "navyblue", "maroon2",
      "lightsalmon1", "lawngreen", "honeydew3", "deeppink4", "cornsilk3", "azure3", "cadetblue2", "gold", "seagreen", "wheat4", "snow3", "purple3",
      "orchid4", "mediumslateblue", "lightsteelblue3", "lightcyan3", "ivory2", "gold2", "darkorange3", "cadetblue4", "azure1", "darkorange1", "paleturquoise2", "steelblue2",
      "tomato1", "seagreen2", "palevioletred4", "navy", "maroon1", "lightsalmon", "lavenderblush4", "honeydew2", "deeppink3", "cornsilk1", "azure2", "cadetblue3",
      "gold4", "seagreen1", "yellow1", "snow1", "purple1", "orchid2", "mediumseagreen", "lightsteelblue2", "lightcyan2", "ivory1", "gold1"
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
#' @param binSize size of bins (bp). Recommend a value bwteen 200 and 10000
#' @param ncores number of cores for parallel computing.
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
ParseReadCounts <- function(data, metaFile = "sample_meta.txt", by = 0.05, prefix = "test", binSize = 1000, ncores = 2) {
  ######################################################
  parParseRaw <- function(x, SEQ) {
    x <- x[x > 0 ]
    totalCPM <- sum(x)
    res <- NULL
    pct <- NULL
    for (K in SEQ) {
      pct <- sum(x[x <= K]) / totalCPM
      res <- c(res, pct)
      if (pct > 0.999) { # exit loop when 99.9% reads have been used
        res <- c(res, rep(1, length(SEQ) - length(res)))
        break
      }
    }
    res
  }
  ######################################################
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
  # Calculate the number of cores
  avai_cores <- detectCores() - 1
  ncores <- ifelse(ncores > avai_cores, avai_cores, ncores)
  if (ncores > 1) {
    cat(paste0("\n\tEnabling parallel computing ( ", ncores, " cores)...\n"))
  }
  tryCatch({
    # Initiate cluster
    cl <- makeCluster(ncores)
    clusterExport(cl, varlist = c("data"), envir = environment())
    # calculate CPM
    CPM <- parLapply(
      cl, 1:ncol(data),
      function(x) data[, x] / sum(data[, x]) * 1000000 * (1000 / binSize)
    )

    CPM <- Reduce("cbind", CPM)
    colnames(CPM) <- colnames(data)
    #  print(dim(CPM))
    #  cat("\n\tCPM was calculateded.\n")

    # remove extreme values potential from centromere regions
    MAX_CPMW <- 150
    CPM[rowSums(CPM > MAX_CPMW) > 0, ] <- 0
    clusterExport(cl, varlist = "CPM", envir = environment())

    # calculate data range and times of loop
    MAX <- 0
    MAX <- parLapply(
      cl, 1:ncol(CPM),
      function(x) max(CPM[, x])
    )

    MAX <- Reduce("cbind", MAX)
    MAX <- as.numeric(ceiling(max(MAX)))
    if (!is.finite(MAX)) {
      MAX <- MAX_CPMW
    }
    SEQ <- seq(0, MAX, by = by)
    # cat("\n\tMAX = ",MAX,"SEQ = ",length(SEQ),"\n")
    if (by > MAX || length(SEQ) < 100) {
      by <- MAX / 100
      SEQ <- seq(0, MAX, by = by)
    }

    clusterExport(cl, varlist = c("SEQ", "CPM", "parParseRaw"), envir = environment())
    options(scipen = 999)
    res <- parLapply(
      cl, 1:ncol(CPM),
      function(x) parParseRaw(CPM[, x], SEQ)
    )
    stopCluster(cl)
    res <- Reduce("cbind", res)
    # cat("\nCPM was parsed.\n")
    colnames(res) <- colnames(CPM)

    dat <- data.frame(cutoff = SEQ, res, check.names = F)
    kept <- rowSums(res) != ncol(res) # delete rows with all 1
    dat <- dat[kept, ]
    output <- paste0(prefix, "_parsedMatrix.txt")
    write.table(dat, output, sep = "\t", quote = F, row.names = F, col.names = T)
    cat("\n\t", output, "[saved]")
  },
  error = function(e) print(e)
  )
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
CalculateSF <- function(data, metaFile = "sample_meta.txt", prefix = "test", xMAX = NA) {
  # calculate scaling factors
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
  data <- data[, c(colnames(data)[1], kept)] # reorder and subset

  MAX_CPM <- max(data[, 1])
  if (!is.na(xMAX)) {
    MAX_CPM <- ifelse(xMAX < MAX_CPM, xMAX, MAX_CPM)
  }
  #---------FUNCTION for quality control-------------------------------------
  QC <- function(data) {
    findLastTurnByPeak <- function(x) {
      # define the last turning point as the highest peak in the density distribution.
      # x must be a vector with names
      # por : proportion of reads
      x <- x[x > 0 & x < 0.99]
      d <- density(x)
      dValues <- diff(d$y)
      turns <- which(dValues[-1] * dValues[-length(dValues)] < 0) + 1
      if (length(turns) == 0) {
        ind <- which.max(d$y)
      } else {
        ind <- max(turns)
      }
      list(por = d$x[ind], cpmw = as.numeric(names(x[x >= d$x[ind]][1])))
    }
    findLastTurnByCutoff <- function(x, cutoff = 0.001) {
      # find the last turning point where dValue is close to cutoff
      # x must be a vector with names
      # por : proportion of reads
      x <- x[x > 0 & x < 0.99]
      dValues <- diff(x)
      dValues <- dValues[dValues > 0]
      cutoff <- ifelse(min(dValues) > cutoff, min(dValues), cutoff)
      dValuesAfterPeak <- dValues[which.max(dValues):length(dValues)]
      ind <- which(dValuesAfterPeak > 0 & dValuesAfterPeak <= cutoff)[1]
      list(por = x[names(dValues)[ind]], cpmw = as.numeric(names(dValues)[ind]))
    }
    findLastTurn <- function(x) {
      # determine the optimal last turning point from two methods
      x <- x[x > 0 & x < 0.99]
      res1 <- findLastTurnByPeak(x)
      res2 <- findLastTurnByCutoff(x)
      delta <- max(x) - res1$por
      cat("\n delta=", delta)
      # to deal with long tail
      if (delta < 0.01) {
        return(res2)
      }
      return(res1)
    }
    findFirstTurnByPeak <- function(x) {
      # define the first turning point before the highest peak
      # x must be a vector with names
      # por : proportion of reads
      x <- x[x > 0 & x < 0.99]
      d1 <- density(x)
      peakX <- d1$x[which.max(d1$y)]
      d2 <- density(d1$x[d1$x < peakX])
      por <- d2$x[which.max(d2$y)]
      list(
        por = por,
        cpmw = as.numeric(names(x[x >= por][1]))
      )
    }

    QC.list <- apply(data[, 2:ncol(data)], 2, FUN = function(x) {
      cpmwCutoff <- 1.20 # empirically define cpmw of input
      names(x) <- data[, 1]
      # find last turning point
      turnLast <- findLastTurn(x)
      turnFirst <- findFirstTurnByPeak(x)
      if (turnLast$cpmw >= cpmwCutoff) {
        QCstr <- "pass"
      } else {
        QCstr <- "fail: complete loss, input or poor enrichment"
      }
      slope <- (turnLast$por - turnFirst$por) / (turnLast$cpmw - turnFirst$cpmw)
      list(
        QC = QCstr, TURNS = paste0(turnFirst$cpmw, ",", turnFirst$por, ",", turnLast$cpmw, ",", turnLast$por),
        xMin = turnFirst$cpmw, yMin = turnFirst$por, xMax = turnLast$cpmw, yMax = turnLast$por, SLOPE = slope
      )
    })
    QC.df <- do.call(rbind.data.frame, QC.list)
    QC.df
  }

  #---------calculate SF-------------------------------------
  QC.df <- QC(data)
  meta <- meta[, !colnames(meta) %in% c("TURNS", "QC", "SF")]
  meta <- cbind(meta, QC.df[rownames(meta), ])
  meta$SF <- NA
  for (ab in unique(meta$ANTIBODY)) {
    inds <- grep(paste0("^", ab, "$"), meta$ANTIBODY) # grep may cause problem sometime
    SF <- round(max(meta$SLOPE[inds]) / meta$SLOPE[inds], 2)
    meta$SF[inds] <- SF
  }
  meta$SF[meta$QC != "pass"] <- NA

  #--------plot each antibody individually------------------------------
  x <- data[, 1]
  imgOutput <- paste0(prefix, "_distribution.pdf")
  pdf(imgOutput, width = 14, height = 8)
  slopesByAb <- NULL
  for (ab in c("All Antibodies", unique(meta$ANTIBODY))) {
    if (ab == "All Antibodies") {
      metaByAb <- meta
    } else {
      metaByAb <- meta[meta$ANTIBODY == ab, ]
    }
    if (length(unique(meta$ANTIBODY)) == 1) {
      # avoid redundant and identical plot
      if (ab == unique(meta$ANTIBODY)) {
        next
      } else {
        ab <- unique(meta$ANTIBODY)
      }
    }
    subsetByAb <- as.data.frame(data[, metaByAb$ID], check.names = F)
    colnames(subsetByAb) <- metaByAb$ID
    # delete rows where all values are out of dataRange 9
    if (ncol(subsetByAb) > 1) {
      kept <- rowSums(subsetByAb > 0.99) != ncol(subsetByAb)
      x <- data[kept, 1]
      subsetByAb <- subsetByAb[kept, ]
    } else {
      x <- data[, 1]
    }
    # Set plot layout
    par(mfrow = c(1, 3), oma = c(0, 0, 5, 0))
    layout.matrix <- matrix(c(1, 2, 3), nrow = 1, ncol = 3)
    layout(
      mat = layout.matrix,
      widths = c(4, 3, 1)
    ) # Widths of the three columns
    par(mar = c(10, 6, 6, 3))
    #-------------plot1: curves--------------------------
    MAX_CPM <- max(x)
    if (ncol(subsetByAb) == 1) {
      totalPages <- 1
    } else {
      totalPages <- 1:ncol(subsetByAb)
    }
    for (r in totalPages) {
      y <- subsetByAb[, r]
      id <- colnames(subsetByAb)[r]
      used <- data.frame(x = x, y = y)
      used <- na.omit(used)
      if (r == 1) {
        plot(used,
          main = "Cumulative Distribution", col = metaByAb$COLOR[r], lwd = 2, xlab = "Cutoff (CPMW)",
          cex.lab = 1.5, cex.axis = 1.5, type = "l", ylab = "Proportion of reads", xlim = c(0, MAX_CPM),
          ylim = c(0.0, 1.1), cex.main = 1.5
        )
      } else {
        lines(used, col = metaByAb$COLOR[r], lty = 1, lwd = 2, pch = 20, cex = 0.1)
      }
      lines(x = c(metaByAb$xMin[r], metaByAb$xMax[r]), y = c(metaByAb$yMin[r], metaByAb$yMax[r]), col = metaByAb$COLOR[r], lty = 3)
    }
    legFontSize <- ifelse(ncol(subsetByAb) < 10, 1.5, 1 + 5 / ncol(subsetByAb))
    legend("bottomright", legend = paste(gsub(".bam", "", metaByAb$ID), paste(", SF=", metaByAb$SF, sep = ""), sep = ""), col = metaByAb$COLOR, pch = 15, bty = "n", ncol = 1, cex = legFontSize)

    #-----------plot2: barplot----------------------------
    cutoffLow <- 3
    cutoffHigh <- 12
    Low <- which.min(abs(x - cutoffLow))
    High <- which.min(abs(x - cutoffHigh))
    barplotDF <- rbind(subsetByAb[Low, ], subsetByAb[High, ] - subsetByAb[Low, ])
    barplotDF <- rbind(barplotDF, 1 - subsetByAb[High, ])
    barplotDF <- as.data.frame(barplotDF, check.names = F)
    colnames(barplotDF) <- metaByAb$ID
    rownames(barplotDF) <- c(
      paste0("< ", cutoffLow),
      paste0(cutoffLow, "~", cutoffHigh),
      paste0("> ", cutoffHigh)
    )

    myCols <- c("grey", "pink", "red")
    xLabels <- colnames(barplotDF)
    legLbls <- rownames(barplotDF)
    par(mar = c(16, 6, 6, 0)) # c(bottom, left, top, right)
    barNum <- ncol(subsetByAb)
    xlimMax <- 1
    barWidth <- ifelse(barNum < 10, xlimMax / barNum * 0.6, xlimMax / barNum * 0.8)
    barSpace <- ifelse(barNum < 10, xlimMax / barNum * 0.4, xlimMax / barNum * 0.2)
    bp <- barplot(as.matrix(barplotDF),
      col = myCols, beside = F, main = "Group by CPMW Range", ylab = "Proportion Of Reads",
      axes = FALSE, axisnames = FALSE, cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5, cex.names = 1.5,
      xlim = c(0, xlimMax), width = barWidth, space = barSpace
    )
    axis(2, cex = 2)
    axis(1, at = bp, labels = FALSE, tck = -0.02)
    xLabels <- gsub(".bam", "", xLabels)
    if (sum(nchar(xLabels) > 20) > 0 || barNum > 10) {
      fontSize <- 1
    } else {
      fontSize <- 1.5
    }

    xLabels <- strtrim(xLabels, 20)
    text(
      x = bp, y = par("usr")[3] - (par("usr")[4] - par("usr")[3]) / 30, labels = xLabels,
      col = metaByAb$COLOR, srt = 60, adj = 1, xpd = TRUE, cex = fontSize
    )

    #---------plot3: legend---------------------------------
    # c(bottom, left, top, right)
    par(mar = c(6, 3, 6, 1), xpd = T)
    plot.new()
    legend("top", fill = rev(myCols), legend = rev(legLbls), ncol = 1, cex = 2, bty = "n", title = "CPMW")
    mtext(ab, outer = TRUE, cex = 1.5, line = 2)
    mtext(imgOutput, outer = TRUE, cex = 1, line = 0)
  }
  garbage <- dev.off()
  cat("\n\t", imgOutput, "[saved]")
  #-----------------------------------------------------
  output <- paste0(prefix, "_SF.txt")
  outDF <- meta[, !colnames(meta) %in% c("SLOPE", "xMin", "yMin", "xMax", "yMax")]
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
    input <- read.table(input, sep = "\t", header = TRUE, fill = TRUE, quote = "", stringsAsFactors = FALSE, check.names = F)
    rownames(input) <- input$ID
  }
  if (!"SF" %in% colnames(input) || !"ANTIBODY" %in% colnames(input) || !"GROUP" %in% colnames(input) || !"QC" %in% colnames(input)) {
    stop("Input looks invalid for BoxplotSF().\n")
  }
  input <- na.omit(input) # exclude samples with SF==NA
  # =======================================================
  plotByAb <- function(metaByGAb, myTitle, PAGE) {
    # no return value, just used side effect to generate plots
    groupSize <- length(unique(metaByAb$GROUP2))
    fontSize <- ifelse(groupSize > 10, 0.6, 1.0)
    if (PAGE == 1) {
      nCol <- ifelse(groupSizeAll > 10, 2, 1)
      marB <- 10
      marL <- 4
    } else {
      nCol <- 1
      marB <- 15
      marL <- 12
    }
    groupLabels <- sort(unique(metaByAb$GROUP2))
    nameOrder <- ordered(metaByAb$GROUP2, levels = groupLabels)
    if (!"COLOR" %in% colnames(metaByAb)) {
      tim10equal <- c("skyblue", "#EF0000", "grey", "#00DFFF", "#50FFAF", "#BFFF40", "#FFCF00", "#FF6000", "#0000FF", "#800000")
      myCols <- tim10equal[as.factor(metaByAb$GROUP2)]
    } else {
      myCols <- metaByAb[!duplicated(metaByAb$GROUP2), "COLOR"]
    }
    #---------panel left---------------
    par(mar = c(marB, marL, marB / 2, 1))
    bp <- boxplot(SF ~ nameOrder,
      data = metaByAb,
      main = myTitle, ylim = c(0, max(metaByAb$SF)),
      las = 2, ylab = "Scaling Factors", xlab = "", xaxt = "n", cex.axis = fontSize, col = myCols,
      cex.lab = 1, cex.main = 1, outline = T, pars = list(outcol = "white")
    )
    axis(1, at = seq_along(groupLabels), labels = FALSE, tck = -0.02)
    text(
      x = seq_along(groupLabels), y = par("usr")[3] - (par("usr")[4] - par("usr")[3]) / 30, srt = 45, adj = 1,
      labels = groupLabels, col = myCols, xpd = TRUE, cex = fontSize
    )
    stripchart(SF ~ nameOrder,
      data = metaByAb,
      vertical = T,
      method = "jitter", add = TRUE, jitter = 0.2, cex = 0.8, pch = 1, col = "#595959"
    )
    #---------panel right---------------

    par(mar = c(marB, 0, marB / 2, 3), xpd = T)
    # c(bottom, left, top, right)
    plot.new()
    legend("top", legend = groupLabels, col = myCols, pch = 15, bty = "n", ncol = nCol, cex = fontSize)
    # mtext(myTitle, outer = TRUE, cex = 1, line = 0)
  }
  # =======================================================

  input$SF <- as.numeric(input$SF)
  input$GROUP2 <- paste(input$GROUP, input$ANTIBODY, sep = ".")
  input <- input[order(input$ANTIBODY), ]
  groupSizeAll <- length(unique(input$GROUP2))
  antibodyList <- unique(input$ANTIBODY)
  imgWidth <- ifelse(groupSizeAll > 15, 8 + groupSizeAll / 10,
    ifelse(length(antibodyList) > 1, 12, 8)
  )
  output <- paste0(prefix, "_boxplot.pdf")
  pdf(output, width = imgWidth, height = 7)
  # page1: All antibodies
  layout.matrix <- matrix(c(1, 2), nrow = 1, ncol = 2)
  layout(layout.matrix, widths = c(4, 3))
  myTitle <- "ScalingFactor ~ All antibodies"
  par(oma = c(2, 2, 2, 2))
  metaByAb <- input
  plotByAb(metaByAb, myTitle, 1)
  if (length(antibodyList) > 1) {
    # avoid redundant and identical plot
    # page2~END: 2 antibodies per page
    num <- 2
    layout.matrix <- matrix(1:(num * 2), nrow = 1, ncol = num * 2)
    layout(layout.matrix, widths = rep(c(4, 3), num))
    for (ab in antibodyList) {
      metaByAb <- input[input$ANTIBODY == ab, ]
      myTitle <- paste0("ScalingFactor ~ ", ab)
      plotByAb(metaByAb, myTitle, 2)
    }
  }
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
#' @param binSize size of bins (bp). Recommend a value bwteen 200 and 10000
#' @param ncores number of cores for parallel computing.
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
ChIPseqSpikeInFree <- function(bamFiles, chromFile = "hg19",
                               metaFile = "sample_meta.txt",
                               prefix = "test", binSize = 1000,
                               ncores = 2,
                               xMAX = NA) {
  # perform ChIP-seq spike-free normalization in one step
  if (binSize < 100 && binSize > 10000) {
    cat(paste0("\n**recommended binSize range 200~10000 (bp); your binSize is", binSize, " **\n"))
  }
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
    rawCountDF <- CountRawReads(bamFiles = bamFiles, chromFile = chromFile, prefix = prefix, binSize = binSize)
  }
  cat("\n\t[--done--]\n")
  cat("\nstep3. parsing raw counts...")
  output2 <- paste0(prefix, "_parsedMatrix.txt")

  # in case of downstream errors, just load existing outFile to avoid re-processing raw counts
  if (file.exists(output2) && disableOverwritten) {
    cat("\n\t", output2, "[just loading the existing file; delete this file first to re-parse rawCount table]")
    parsedDF <- read.table(output2, sep = "\t", header = TRUE, fill = TRUE, stringsAsFactors = FALSE, quote = "", row.names = NULL, check.names = F)
  } else {
    steps <- 0.05 * binSize / 1000
    parsedDF <- ParseReadCounts(data = rawCountDF, metaFile = meta, prefix = prefix, by = steps, binSize = binSize, ncores = ncores)
  }
  cat("\n\t[--done--]\n")
  cat("\nstep4. calculating scaling factors...")
  result <- CalculateSF(data = parsedDF, metaFile = meta, prefix = prefix, xMAX = xMAX)
  cat("\t[--done--]\n")
  cat("\nstep5. ploting scaling factors...")
  BoxplotSF(result, prefix)
  cat("\n\t[--done--]\n")
  invisible(result)
}
