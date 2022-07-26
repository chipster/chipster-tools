# TOOL dada2-filter.R: "Filter input sequences with dada2" (Given a tar package of fastg files, this tool filters the input sequences which don't fullfill the user defined criteria. For more information see the manual)
# INPUT reads.tar: "Tar package containing the FASTQ files" TYPE GENERIC
# OUTPUT OPTIONAL summary.tsv
# OUTPUT OPTIONAL summary.txt
# OUTPUT filtered.tar
# PARAMETER OPTIONAL truncf: "Truncate forward reads after this amount of bases" TYPE INTEGER FROM 0 DEFAULT 0 (Default 0 means no truncation. Truncate reads after truncLen bases. Reads shorter than this are discarded.) 
# PARAMETER OPTIONAL truncr: "Truncate reverse reads after this amount of bases" TYPE INTEGER FROM 0 DEFAULT 0 (Default 0 means no truncation. Truncate reads after truncLen bases. Reads shorter than this are discarded.) 
# PARAMETER OPTIONAL maxns: "Discard input sequences with more than specified number of Ns" TYPE INTEGER FROM 0 DEFAULT 0 (Sequences with more than specified number of Ns will be discarded. Note that dada does not allow any Ns.)
# PARAMETER OPTIONAL maxeef: "Discard forward sequences with more than the specified number of expected errors" TYPE DECIMAL FROM 0 (After truncation, reads with higher than maxEE "expected errors" will be discarded. If this parameter is not set, no expected error filtering is done.)
# PARAMETER OPTIONAL maxeer: "Discard reverse sequences with more than the specified number of expected errors" TYPE DECIMAL FROM 0 (After truncation, reads with higher than maxEE "expected errors" will be discarded. If this parameter is not set, no expected error filtering is done.)
# PARAMETER OPTIONAL truncq: "Truncuate sequence after this base quality" TYPE INTEGER FROM 0 (Truncate reads at the first instance of a quality score less than or equal to the specified number.)
# RUNTIME R-4.1.1

# ES 15.07.2022
source(file.path(chipster.common.path,"tool-utils.R"))
source(file.path(chipster.common.path,"zip-utils.R"))

# load library dada2
library(dada2)
#packageVersion("dada2")

#check out if the file is compressed and if so unzip it
unzipIfGZipFile("reads.tar")

# Read the contents of the tar file into a list
system("tar tf reads.tar > tar.contents")
file.list <- scan("tar.contents",what = "",sep = "\n")

# Check that the input is a valid tar file
if (length(file.list) == 0) {
  stop(paste('CHIPSTER-NOTE: ',"It seems your input file is not a valid Tar package. Please check your input file."))
}
# Check if tar packa contains folders
if (grepl("/",file.list[1])) {
  stop(paste('CHIPSTER-NOTE: ',"It seems your Tar package contains folders. The FASTQ files need to be in the root of the package, not in subfolders."))
}

# Make input and output folders 
system("mkdir input_folder")
system("mkdir output_folder")

# untar the tar package to input_folder and list the filenames
untar("reads.tar", exdir = "input_folder")
filenames <- list.files("input_folder", full.names=TRUE)

# Sort the filenames, samples have different names
filenames <- sort(filenames)

# check if the lenght of files in the input_folder is even, else error
number <- length(filenames)%%2
if (number != 0 && length(filenames)<2){
    stop(paste('CHIPSTER-NOTE: ',"It seems that some of your fastq files doesn`t have a pair"))
    }

# put the forward files to fnFs, assume that forward reads have the same name but different tag than reverse reads
# forward reads should be before reverse reads after sort() function
forward <- seq(1,length(filenames),by=2)
fnFs <- filenames[forward]

reverse <- seq(2,length(filenames), by=2)
fnRs <- filenames[reverse]

# put the reverse and forward reads to own folders  problem with pattern
#fnFs <- sort(list.files("input_folder", pattern="_R1_001.fastq", full.names = TRUE))
#fnRs <- sort(list.files("input_folder", pattern="_R2_001.fastq", full.names = TRUE))

# take out the sample names splitting with _
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)


# make a folder for the reverse and forward reads for the filterandtrim function. Rename output files with F=forward and R=reverse
filtFs <- file.path("output_folder", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path("output_folder", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# if maxeer or maxeef not selected 
if (is.na(maxeef)){
  maxeef='inf'
}
if (is.na(maxeer)){
  maxeer='inf'
}
# run command filterAndTrim
# if truncq not selected
if (is.na(truncq)){
  out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(truncf,truncr),
              maxN=maxns, maxEE=c(maxeef,maxeer), rm.phix=TRUE,
              compress=TRUE, multithread=TRUE, verbose=TRUE)
             
  }else{
  out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(truncf,truncr),
              maxN=maxns, maxEE=c(maxeef,maxeer), truncQ=truncq, rm.phix=TRUE,
              compress=TRUE, multithread=TRUE, verbose=TRUE)}

sink("summary.txt")
	cat("\n\n\n")
  print(out)
	cat("\n\n\n")
sink()

# make summary.tsv table
rownames(out) <- sample.names
write.table(out, file ="summary.tsv", sep='\t')

#make a output tar package named filtered.tar and gzip
system("gzip output_folder/*.fq")
system("cd output_folder && tar cf ../filtered.tar *")

#EOF