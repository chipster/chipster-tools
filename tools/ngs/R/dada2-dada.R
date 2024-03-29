# TOOL dada2-dada.R: "Sample inference" (Given a tar package of FASTQ files, this tool runs the learnErrors and dada commands from the DADA2 library. If the visualize error rates parameter is set to yes, the error rates are visualized to a pdf file. For the dada function the ambigious bases Ns needs to be removed before from the fastq files. This tool can be used either for single or paired end reads. Please check the manual for more information.)
# INPUT reads.tar: "Tar package containing the FASTQ files" TYPE GENERIC
# OUTPUT dada_forward.Rda
# OUTPUT OPTIONAL dada_reverse.Rda
# OUTPUT summary.txt
# OUTPUT OPTIONAL plotErrors.pdf
# PARAMETER type: "Is the data paired end reads, IonTorrent or other types of single end reads?" TYPE [paired, IonTorrent, single] DEFAULT paired (If the reads were produced with Ion Torrent, DADA2 setting for that techonology will be applied.)
# PARAMETER ploterr: "Do you want to visualize the estimated error rates?" TYPE [yes,no] DEFAULT no (Do you want to visualize the error rates to a pdf file)
# PARAMETER pool: "Type of pooling" TYPE [independent, pseudo-pooling] DEFAULT independent (If this is set to pseudo-pooling, the dada algorithm will perform independent processign twice, which makes the sensitivity better but processing takes twice longer. Check the manual.)
# RUNTIME R-4.1.1-asv


# ES 15.07.2022
# SLOTS 2  #with slots2 tool is faster ~2x
# OUTPUT OPTIONAL log.txt
# OUTPUT OPTIONAL log2.txt
# OUTPUT OPTIONAL summary_stats.tsv

source(file.path(chipster.common.lib.path, "tool-utils.R"))
source(file.path(chipster.common.lib.path, "zip-utils.R"))

# load library dada2
library(dada2)
# packageVersion("dada2")
library(ggplot2) # just for renaming plot title
# Read the contents of the tar file into a list
system("tar tf reads.tar > tar.contents")
file.list <- scan("tar.contents", what = "", sep = "\n")

# Check that the input is a valid tar file
if (length(file.list) == 0) {
  stop(paste("CHIPSTER-NOTE: ", "It seems your input file is not a valid Tar package. Please check your input file."))
}
# Check if tar packa contains folders
if (grepl("/", file.list[1])) {
  stop(paste("CHIPSTER-NOTE: ", "It seems your Tar package contains folders. The FASTQ files need to be in the root of the package, not in subfolders."))
}

# Make input folder
system("mkdir input_folder")

# untar the tar package to input_folder and list the filenames
untar("reads.tar", exdir = "input_folder")
filenames <- list.files("input_folder", full.names = TRUE)

# Sort the filenames _R1 _R2, samples have different names
filenames <- sort(filenames)

# if the reads are paired end
if (type == "paired") {
  # check if the lenght of files in the input_folder is even (all the fastq files have a pair), else error
  number <- length(filenames) %% 2
  if (number != 0) {
    stop(paste("CHIPSTER-NOTE: ", "It seems that some of your fastq files doesn`t have a pair"))
  }
  # put the forward files to fnFs, assume that forward reads have the same name but different tag than reverse reads
  # forward reads should be before reverse reads after sort() function, should use filterAndTrim() before, if there is a problem.
  forward <- seq(1, length(filenames), by = 2)
  fnFs <- filenames[forward]

  reverse <- seq(2, length(filenames), by = 2)
  fnRs <- filenames[reverse]

  # put the reverse and forward reads to own folders  problem with pattern
  # fnFs <- sort(list.files("input_folder", pattern="_R1_001.fastq", full.names = TRUE))
  # fnRs <- sort(list.files("input_folder", pattern="_R2_001.fastq", full.names = TRUE))

  # take out the sample names
  sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

  sink(file = "log.txt")
  # Run learnErrors (slow)
  set.seed(100)
  cat("\nLearn error rates of forward reads:\n")
  errF <- learnErrors(fnFs, multithread = as.integer(chipster.threads.max))
  cat("\nLearn error rates of reverse reads:\n")
  errR <- learnErrors(fnRs, multithread = as.integer(chipster.threads.max))
  sink()

  # make a pdf file for visualizing error rates if ploterrors parameter is yes
  if (ploterr == "yes") {
    pdf("plotErrors.pdf", , width = 13, height = 7)
    print(plotErrors(errF, nominalQ = TRUE) + labs(title = "Estimated error rates of forward reads"))
    print(plotErrors(errR, nominalQ = TRUE) + labs(title = "Estimated error rates of reverse reads"))
    dev.off()
  }

  # run dereplicate, no need dada() can handle fastq files, same result
  # derepF1 <- derepFastq(fnFs)
  # derepR1 <- derepFastq(fnRs)

  # run dada() and make a summary. If pool=pseudo then run with parameter pseudo
  if (pool == "independent") {
    sink(file = "log2.txt")
    cat("\nDada function of forward reads:\n")
    dadaFs <- dada(fnFs, err = errF, multithread = as.integer(chipster.threads.max))
    cat("\nDada function of reverse reads:\n")
    dadaRs <- dada(fnRs, err = errR, multithread = as.integer(chipster.threads.max))
    cat("\n")
    cat("\nDada-class objects describing DADA2 denoising results of forward reads:\n\n")
    print(dadaFs)
    cat("\nDada-class objects describing DADA2 denoising results of reverse reads:\n\n")
    print(dadaRs)
    sink()
  } else {
    sink(file = "log2.txt")
    cat("\nDada function of forward reads:\n")
    dadaFs <- dada(fnFs, err = errF, pool = "pseudo", multithread = as.integer(chipster.threads.max))
    cat("\nDada function of reverse reads:\n")
    dadaRs <- dada(fnRs, err = errR, pool = "pseudo", multithread = as.integer(chipster.threads.max))
    cat("\n")
    cat("\nDada-class objects describing DADA2 denoising results of forward reads:\n\n")
    print(dadaFs)
    cat("\nDada-class objects describing DADA2 denoising results of reverse reads:\n\n")
    print(dadaRs)
    sink()
  }
  # save the dada-class objects to .Rda file
  save(dadaFs, file = "dada_forward.Rda")
  save(dadaRs, file = "dada_reverse.Rda")

  # ----------------- if single end reads, the same thing without reverse parts, if ionTorrent change options --------------------------------------------------
} else {
  # take out the sample names
  sample.names <- sapply(strsplit(basename(filenames), "_"), `[`, 1)

  sink(file = "log.txt")
  # Run learnErrors (slow)
  cat("\nLearn error rates:\n")
  set.seed(100)
  errF <- learnErrors(filenames, multithread = as.integer(chipster.threads.max))
  sink()

  # make a pdf file for visualizing error rates if ploterrors parameter is yes
  if (ploterr == "yes") {
    pdf("plotErrors.pdf", , width = 13, height = 7)
    print(plotErrors(errF, nominalQ = TRUE) + labs(title = "Estimated error rates"))
    dev.off()
  }

  # run dada() and make a summary. If pool=pseudo then run with parameter pseudo
  if (pool == "independent") {
    sink(file = "log2.txt")
    cat("\n# Dada() function\nType of pooling: Independent\n")
    if (type == "IonTorrent") { # set homopolymer gab penalty -1 and band size 32
      setDadaOpt(HOMOPOLYMER_GAP_PENALTY = -1, BAND_SIZE = 32)
      cat("Data: IonTorrent\nHOMOPOLYMER_GAB_PENALTY= ", getDadaOpt("HOMOPOLYMER_GAP_PENALTY"), "\n")
      dadaFs <- dada(filenames, err = errF, multithread = as.integer(chipster.threads.max))
    } else {
      cat("Data: single end reads\n")
      dadaFs <- dada(filenames, err = errF, multithread = as.integer(chipster.threads.max), verbose = TRUE)
    }
    cat("\n")
    cat("\nDada-class objects describing DADA2 denoising results:\n\n")
    print(dadaFs)
    sink()
  } else {
    sink(file = "log2.txt")
    cat("\n# Dada() function\nType of pooling: Pseudo-pooling\n")
    if (type == "IonTorrent") {
      setDadaOpt(HOMOPOLYMER_GAP_PENALTY = -1, BAND_SIZE = 32)
      cat("Data: IonTorrent\nHOMOPOLYMER_GAB_PENALTY= ", getDadaOpt("HOMOPOLYMER_GAP_PENALTY"), "\n")
      dadaFs <- dada(filenames, err = errF, multithread = as.integer(chipster.threads.max))
    } else {
      cat("Data: single end reads\n")
      dadaFs <- dada(filenames, err = errF, pool = "pseudo", multithread = as.integer(chipster.threads.max))
    }
    cat("\n")
    cat("\nDada-class objects describing DADA2 denoising results:\n\n")
    print(dadaFs)
    sink()
  }
  # save the dada-class objects to .Rda file
  save(dadaFs, file = "dada_forward.Rda")
}

# merge the log files to one file, remove few rows
file.create("summary.txt")
rows <- readLines("log.txt")
for (row in rows) {
  write(row, file = "summary.txt", append = TRUE)
}
rows <- readLines("log2.txt")
for (row in rows) {
  if (grepl("Key parameters:", row)) {
    write(paste0("\nParameters used for the DADA() function:\n", row), file = "summary.txt", append = TRUE)
    break
  }
}
for (row in rows) {
  if (!grepl("dada-class:", row) && !grepl("Key ", row) && !grepl("   selfConsist", row)) {
    write(row, file = "summary.txt", append = TRUE)
  }
}

# EOF
