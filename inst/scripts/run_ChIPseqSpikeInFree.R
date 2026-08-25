#!/usr/bin/env Rscript
# last update: 8/25/2026 by Hongjian Jin @ St Jude Children's Research Hospital
######################################################
Rver <- as.numeric(paste0(R.version$major, ".", floor(as.numeric(R.version$minor))))
if (Rver < 3.6) {
  cat("\n**please used R/3.6.1 module to run ChIPseqSpikeInFree**\n")
  stop(paste0("Your current R version is ", paste0(R.version$major, ".", R.version$minor), "\n\n"))
}

######################################################
ChIPseqSpikeInFree2 <- function(bamFiles, chromFile = genome,
                               metaFile = metaFile,
                               prefix = prefix, binSize = binSize,
                               cutoff_QC = cutoffQC, maxLastTurn=fracMax,
                               minFirstTurn=fracMin, ncores =ncores, singleEnd=FALSE
                               ) 
 {
   # perform ChIP-seq spike-free normalization in one step
  if (binSize < 200 || binSize > 10000) {
    stop(paste0("\n**recommended binSize range 200~10000 (bp); your binSize is", binSize, " **\n"))
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
    rawCountDF <- CountRawReads(bamFiles = bamFiles, chromFile = chromFile,  singleEnd = singleEnd, prefix = prefix, binSize = binSize)
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
  SF <- CalculateSF(data = parsedDF, metaFile = meta, minFirstTurn = minFirstTurn, maxLastTurn=maxLastTurn, cutoff_QC=cutoff_QC)
  cat("\nstep5. ploting distribution curves and bars ...")
  PlotDistr(parsedDF, SF, prefix, xlimMaxCPMW=xlimMaxCPMW )
  cat("\t[--done--]\n")
  cat("\nstep6. ploting scaling factors...")
  BoxplotSF(SF, prefix)
  cat("\n\t[--done--]\n")
  invisible(SF)
}

######################################################

suppressMessages(library(optparse))
options(warn=-1)
option_list <- list(
  make_option(c("-L", "--bamList"),
    type = "character", default = NA,
    help = "bam.lst, the same bam file list used for counting raw reads"
  ),
  make_option(c("-G", "--genome"),
    type = "character", default = NA,
    help = "mm9, mm10, hg19, hg38 are supported; or a chrom.size file"
  ),
  make_option(c("-m", "--metaData"),
    type = "character", default = NA,
    help = "meta file with at least three required columns:\n\t\tID,ANTIBODY and GROUP. [COLOR is optional]."
  ), make_option(c("-o", "--prefix"),
    type = "character", default = NA,
    help = "prefix for output filename"
  ),
  make_option(c("-b", "--binSize"),
    type = "numeric", default = 1000,
    help = "bin size for read counting.Recommend a value between 200 and 10000 (bp) [default %default]."
  ),
  make_option(c("-r", "--replot"),
    type = "logical", default = F,
    help = "FALSE - standard mode; TRUE - replotting mode [default %default]."
  ),
    make_option(c("-E", "--paired"),
    type = "logical", default = F,
    help = "FALSE - single-end; TRUE - paired-end [default %default]."
  ),
 make_option(c("-f", "--minFirstTurn"),
    type = "numeric", default = NA,
    help = "numeric; define minFirstTurn [between 0.20 and 0.8], i.e. exclude this fraction of noisy reads before estimating the first turning point [default %default]."
  ),

  make_option(c("-l", "--maxLastTurn"),
    type = "numeric", default = 0.99,
    help = "numeric; define maxLastTurn [between 0.90 and 0.99], i.e. fraction of reads used to estimate the last turning point [default %default]."
  ),
  make_option(c("-c", "--cutoffQC"),
    type = "numeric", default = 1.2,
    help = "numeric; define cutoff for QC [between 0.95 and 1.2]. smaller value, looser cutoff[default %default]."
  ),
    make_option(c("-x", "--xmax"),
    type = "numeric", default = NULL,
    help = "numeric; define xlim_max value in distribution plot[default %default]."
  ),
  make_option(c("-n", "--ncores"),
    type = "numeric", default = 2,
    help = "number of cores for parallel computing [default %default]."
  )
)

parser <- OptionParser(
  usage = paste("\n[--Standard Mode--]",
                "%prog -L <bam.lst> -G <genome or chrom.size> 
                -m <sample_meta.txt> -o <output_prefix> [-f <minFirstTurn> 
                 -l <maxLastTurn> -c <cutoffQC> -b <bin_size> -n <ncores> ]",
                "\n[--Replotting Mode--]",
                " %prog -m <sample_meta.txt> -o <output_prefix> -r TRUE ", sep="\n"),
  description = "\nCalculate ChIPseqSpikeFree normalization factors and optionally replot them",
  option_list = option_list,
  epilogue = paste("Standard mode: Arguments BAMLIST, GENOME ,METADATA, and OUTPUT are required; others are optional.\n",
                    "Replotting mode: Arguments -m, -o and -r TRUE are required; read-counting options -L, -G, -b, -E and -n are ignored.\n", 
                    "\nExample:",
                    " %prog -L bam.lst -G hg19 -m sample_meta.txt -o test",
                    " %prog -L bam.lst -G hg19 -m sample_meta.txt -l 0.97 -c 1 -o test [for H3K9me3]",
                    " %prog -m sample_meta.csf.txt -o H3K9me3.csf -l 0.97 -c 1 -r TRUE ",
                    " %prog -m sample_meta.csf.txt -o H3K9me3.csf -x 10 -r TRUE",
                    " %prog -m sample_meta.csf.txt -o H3K9me3.csf -f 0.5 -l 0.97 -c 1 -r TRUE \n\n",
                    sep = "\n")
)
options(warn=0)

args <- NA
result <- tryCatch({
  args <- parse_args(parser, positional_arguments = TRUE) # TRUE = c(0, Inf), FALSE,1,  c(1,2)
}, warning = function(w) {
  message(w)
  cat("\n\n")
}, error = function(e) {
  message(e)
  cat("\n\n")
})
opt <- args$options
if (is.na(opt$metaData) || is.na(opt$prefix)) {
  print_help(parser)
  quit("no")
}
metaFile <- normalizePath(opt$metaData)
prefix <- opt$prefix
binSize <- opt$binSize
ncores <- opt$ncores
replot <- opt$replot
paired <- opt$paired
if (is.na(opt$minFirstTurn)){
   fracMin <- "auto"
}else{
   fracMin <- opt$minFirstTurn
}
fracMax <- opt$maxLastTurn
cutoffQC <- opt$cutoffQC
xlimMaxCPMW <- opt$xmax
suppressMessages(library(ChIPseqSpikeInFree))

######################################################
if (replot){
    cat("\n[--Replotting Mode--]\n")
    if (is.na(metaFile) || is.na(prefix)){
        stop("\n Arguments METADATA, OUTPUT and REPLOT are required in replotting mode.\n\n")
    }  
    meta <- ReadMeta(metaFile)
    parsedMatrix <- paste0(prefix,"_parsedMatrix.txt")
    if (!file.exists(parsedMatrix)) {
      stop(paste0("parsedMatrix file not found:", parsedMatrix,"\n\n"))
    }
    parsedDF <- read.table(parsedMatrix, sep = "\t", header = TRUE, fill = TRUE, stringsAsFactors = FALSE, quote = "", row.names = NULL, check.names = F)
    SF <- CalculateSF(parsedDF, metaFile = meta, minFirstTurn=fracMin, maxLastTurn=fracMax, cutoff_QC=cutoffQC)
    PlotDistr(parsedDF, SF, paste0(prefix,"_replot"), xlimMaxCPMW=xlimMaxCPMW )
    BoxplotSF(SF, paste0(prefix,"_replot"))
    cat("\nCheers!\n\n")
    quit("no")
}
######################################################
cat("\n[--Standard Mode--]\n")

bamList <- opt$bamList
genome <- opt$genome

if (!file.exists(bamList)) {
  stop(paste0("bamList file not found:", bamList))
}
if (!tolower(genome) %in% c("mm9", "mm10", "hg19", "hg38")) {
  if (!file.exists(genome)) {
    stop(paste0("**Warning: you may use an unsupported genome or \n
    your chrom.size file doesn\'t exist,", genome, "**\n"))
  }
}

bams <- scan(file = bamList, character(), quote = "", quiet = T)
if (paired){
   ChIPseqSpikeInFree2(bamFiles = bams, chromFile = genome, metaFile = metaFile, prefix = prefix, maxLastTurn=fracMax, minFirstTurn=fracMin, cutoff_QC=cutoffQC, ncores=ncores, binSize=binSize)
   # ChIPseqSpikeInFree(bamFiles = bams, chromFile = genome, metaFile = metaFile, prefix = prefix, maxLastTurn=fracMax, cutoff_QC=cutoffQC, ncores=ncores, binSize=binSize)
}else{
   ChIPseqSpikeInFree2(bamFiles = bams, chromFile = genome, metaFile = metaFile, prefix = prefix, maxLastTurn=fracMax, minFirstTurn=fracMin, cutoff_QC=cutoffQC, ncores=ncores, binSize=binSize,singleEnd=TRUE)
}

cat("\nCheers!\n\n")

######################################################
# last updated: 8/25/2026 PM
######################################################
# 1/27/2021, added paired-end mode
# 10/30/2025, exposed -x --xmax





