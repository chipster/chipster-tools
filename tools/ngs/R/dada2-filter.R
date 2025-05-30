# TOOL dada2-filter.R: "Filter and trim sequences with DADA2" (Given a tar package of FASTQ files, this tool filters the input sequences which don't fullfill the user defined criteria. This tool can be used either for single or paired end reads. If the reads are single end, then use only the parameters for forward reads. For more information please check the manual.)
# INPUT reads.tar: "Tar package containing the FASTQ files" TYPE GENERIC
# INPUT OPTIONAL input_list.txt: "List of FASTQ files by sample" TYPE GENERIC (If the FASTQ files are not assigned into samples correctly, you can give a file containing this information. Check instructions from manual)
# OUTPUT filtered.fastqs.tar
# OUTPUT summary.tsv
# OUTPUT samples.fastqs.txt
# PARAMETER paired: "Is the data paired end or single end reads" TYPE [paired, single] DEFAULT paired (If single end reads, use only those forward parameters.)
# PARAMETER OPTIONAL truncf: "Truncate forward reads after this amount of bases" TYPE INTEGER FROM 0 DEFAULT 0 (Default 0 means no truncation. Truncate reads after truncLen bases. Reads shorter than this are discarded. You can use this parameter for single and paired end reads.)
# PARAMETER OPTIONAL truncr: "Truncate reverse reads after this amount of bases" TYPE INTEGER FROM 0 DEFAULT 0 (Default 0 means no truncation. Truncate reads after truncLen bases. Reads shorter than this are discarded. Use only for paired end reads.)
# PARAMETER OPTIONAL maxns: "Discard input sequences with more than specified number of Ns" TYPE INTEGER FROM 0 DEFAULT 0 (Sequences with more than the specified number of Ns will be discarded. Note that the dada function does not allow any Ns.)
# PARAMETER OPTIONAL maxeef: "Discard forward sequences with more than the specified number of expected errors" TYPE DECIMAL FROM 0 (After truncation, reads with more than this amount of expected errors will be discarded. If this parameter is not set, no expected error filtering is done. You can use this parameter for single and paired end reads.)
# PARAMETER OPTIONAL maxeer: "Discard reverse sequences with more than the specified number of expected errors" TYPE DECIMAL FROM 0 (After truncation, reads with more than this amount of expected errors will be discarded. If this parameter is not set, no expected error filtering is done. Use only for paired end reads.)
# PARAMETER OPTIONAL truncq: "Truncate reads after this base quality" TYPE INTEGER FROM 0 DEFAULT 2 (Truncate reads at the first instance of a quality score less than or equal to the specified number. Setting this parameter to 0, turns this behaviour off.)
# PARAMETER OPTIONAL triml: "The number of nucleotides to remove from start of each read" TYPE INTEGER FROM 0 DEFAULT 0 (The number of nucleotides to remove from the start of each read. If both truncLen and trimLeft are provided, filtered reads will have length truncLen - trimLeft.)
# PARAMETER OPTIONAL minlen: "Remove reads which are shorter than this" TYPE INTEGER FROM 0 DEFAULT 0 (Removes reads which are shorter than the specified value. Min length is enforced after all other trimming and truncation. This parameter is especially usefull when truncLen parameter is not used for example with ITS data.)
# RUNTIME R-4.1.1-asv
# SLOTS 2


# ES 15.07.2022
# OUTPUT OPTIONAL summary.txt

source(file.path(chipster.common.lib.path, "tool-utils.R"))
source(file.path(chipster.common.lib.path, "zip-utils.R"))

# load library dada2
library(dada2)
# packageVersion("dada2")

# check out if the file is compressed and if so unzip it
unzipIfGZipFile("reads.tar")

# print(chipster.threads.max)
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

# Make input and output folders
system("mkdir input_folder")
system("mkdir output_folder")

# untar the tar package to input_folder and list the filenames
untar("reads.tar", exdir = "input_folder")

# list the full file names
filenames <- list.files("input_folder", full.names = TRUE)



## ---------------------------------
if (paired == "paired") {
  # Use input_list if provided. Else use tar package file names tool to generate the input list
  if (fileOk("input_list.txt")) {
    txt_filenames <- c() # name of the files in the txt file
    sample.names <- c()
    input <- readLines("input_list.txt")
    # take out the file names and sample name and put them to one vector, those are separeted with '\t'
    for (row in input) {
      sample <- strsplit(row, "\t", fixed = TRUE)
      sample_name <- trimws(sample[[1]][1])
      sample.names <- c(sample.names, sample_name) # sample names first row
      txt_filenames <- c(txt_filenames, trimws(sample[[1]][2])) # forward read
      txt_filenames <- c(txt_filenames, trimws(sample[[1]][3])) # reverse read
    }
    # if everything fine change the filenames variable and use it, add also the folder name: input_folder/
    if (length(txt_filenames) != length(filenames)) {
      stop(paste("CHIPSTER-NOTE: ", "It seems that the list of FASTQ files .txt file has different amount of filenames than the .tar package"))
    } else {
      filenames <- paste0("input_folder/", txt_filenames)
    }
  } else {
    # Sort the filenames from tar package, samples have different names
    filenames <- sort(filenames)
    txt_filenames <- list.files("input_folder")
  }

  # check if the lenght of files in the input_folder is even, else error
  number <- length(filenames) %% 2
  if (number != 0) {
    stop(paste("CHIPSTER-NOTE: ", "It seems that some of your fastq files doesn`t have a pair"))
  }

  # put the forward files to fnFs, assume that forward reads have the same name but different tag than reverse reads
  # forward reads should be before reverse reads after sort() function
  forward <- seq(1, length(filenames), by = 2)
  fnFs <- filenames[forward]

  reverse <- seq(2, length(filenames), by = 2)
  fnRs <- filenames[reverse]


  # take out the sample names splitting with _ if no input_file
  if (!fileOk("input_list.txt")) {
    sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
  }

  # put the reverse and forward reads to own folders  problem with pattern
  # fnFs <- sort(list.files("input_folder", pattern="_R1_001.fastq", full.names = TRUE))
  # fnRs <- sort(list.files("input_folder", pattern="_R2_001.fastq", full.names = TRUE))


  # make a folder for the reverse and forward reads for the filterandtrim function. Rename output files with F=forward and R=reverse
  filtFs <- file.path("output_folder", paste0(sample.names, "_F_filt.fastq.gz"))
  filtRs <- file.path("output_folder", paste0(sample.names, "_R_filt.fastq.gz"))
  names(filtRs) <- sample.names
  names(filtFs) <- sample.names


  # if maxeer or maxeef not selected
  if (is.na(maxeef)) {
    maxeef <- "inf"
  }
  if (is.na(maxeer)) {
    maxeer <- "inf"
  }
  # run filterAndTrim
  out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
    truncLen = c(truncf, truncr),
    maxN = maxns, maxEE = c(maxeef, maxeer), truncQ = truncq, minLen = minlen, rm.phix = TRUE,
    compress = TRUE, multithread = as.integer(chipster.threads.max), verbose = TRUE
  )

  if (length(filenames) > 2) {
    # print(length("double_many",filenames)
    rownames(out) <- sample.names
    write.table(out, file = "summary.tsv", sep = "\t")
  } else {
    # print(length("pair",filenames))
    write.table(out, file = "summary.tsv", sep = "\t", row.names = FALSE)
    print(out)
  }
  # make a summary file
  # sink("summary.txt")
  # 	cat("\n\n\n")
  #  print(out)
  # 	cat("\n\n\n")
  # sink()

  # make a samples.fastqs.txt file to show how fastq files were assigned to samples, every second file forward
  file.create("samples.fastqs.txt")
  x <- 1
  for (name in sample.names) {
    line <- paste(name, "\t", txt_filenames[x], "\t", txt_filenames[x + 1])
    x <- x + 2
    write(line, file = "samples.fastqs.txt", append = TRUE)
  }
  ## ----------------------------------------------------------------------------------
} else { # single reads, otherwise is the same. Has just forward reads.
  # Use input_list if provided. Else use tar package file names tool to generate the input list

  if (fileOk("input_list.txt")) {
    txt_filenames <- c() # name of the files in the txt file
    sample.names <- c()
    input <- readLines("input_list.txt")
    # take out the file names and sample name and put them to one vector, those are separeted with '\t'
    for (row in input) {
      sample <- strsplit(row, "\t", fixed = TRUE)
      sample_name <- trimws(sample[[1]][1])
      sample.names <- c(sample.names, sample_name) # sample names first row
      txt_filenames <- c(txt_filenames, trimws(sample[[1]][2])) # forward read
    }

    # if everything fine change the filenames variable and use it, add also the folder name: input_folder/
    if (length(txt_filenames) != length(filenames)) {
      stop(paste("CHIPSTER-NOTE: ", "It seems that the list of FASTQ files .txt file has different amount of filenames than the .tar package"))
    } else {
      filenames <- paste0("input_folder/", txt_filenames)
    }
  } else {
    # Sort the filenames from tar package, samples have different names, take put sample names _
    filenames <- sort(filenames)
    txt_filenames <- list.files("input_folder")
    sample.names <- sapply(strsplit(basename(filenames), "_"), `[`, 1)
  }
  filtreads <- file.path("output_folder", paste0(sample.names, "_filt.fastq.gz"))
  # if maxeer or maxeef not selected
  if (is.na(maxeef)) {
    maxeef <- "inf"
  }

  # run filterAndTrim
  out <- filterAndTrim(filenames, filtreads,
    truncLen = truncf,
    maxN = maxns, maxEE = maxeef, truncQ = truncq, minLen = minlen, rm.phix = TRUE,
    compress = TRUE, multithread = as.integer(chipster.threads.max), verbose = TRUE
  )

  if (length(filenames) > 1) {
    # print(length("sing_many",filenames))
    rownames(out) <- sample.names
    write.table(out, file = "summary.tsv", sep = "\t")
  } else {
    # print(length(filenames))
    write.table(out, file = "summary.tsv", sep = "\t", row.names = FALSE)
  }
  x <- 1
  file.create("samples.fastqs.txt")
  for (name in sample.names) {
    line <- paste(name, "\t", txt_filenames[x], "\t")
    x <- x + 1
    write(line, file = "samples.fastqs.txt", append = TRUE)
  }
}

# make summary.tsv table
# rownames(out) <- sample.names
# write.table(out, file ="summary.tsv", sep='\t', row.names=FALSE)

# make a output tar package named filtered.tar and gzip
system("gzip output_folder/*.fq")
system("cd output_folder && tar cf ../filtered.fastqs.tar *")

# EOF
