# TOOL dada2-filter.R: "Filter input reads with dada2"
# INPUT reads.tar: "Tar package containing the FASTQ files" TYPE GENERIC
# OUTPUT OPTIONAL summary.tsv
# OUTPUT OPTIONAL summary.txt
# OUTPUT filtered.tar
# PARAMETER OPTIONAL truncf: "Base truncate forward after" TYPE INTEGER FROM 0 DEFAULT 0 (0 means no )
# PARAMETER OPTIONAL truncr: "Base truncate reverse after" TYPE INTEGER FROM 0 DEFAULT 0
# PARAMETER OPTIONAL maxns: "Discard input sequences with more than this N" TYPE INTEGER FROM 0 DEFAULT 0
# PARAMETER OPTIONAL maxeer: "Discard reverse sequences with more than the specified number of expected errors" TYPE DECIMAL FROM 0 
# PARAMETER OPTIONAL maxeef: "Discard forward sequences with more than the specified number of expected errors" TYPE DECIMAL FROM 0 
# PARAMETER OPTIONAL truncq: "Truncuate sequence after this base quality" TYPE INTEGER FROM 0 DEFAULT 2
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
filenames <- list.files("input_folder")

 # Sort the filenames _R1 _R2, samples have different names
filenames <- sort(filenames)

# check if the lenght of files in the input_folder is even (all the fastq files have a pair), else error
number <- length(filenames)%%2
if (number != 0 && length(filenames)<2){
    stop(paste('CHIPSTER-NOTE: ',"It seems that some of your fastq files doesn`t have a pair"))
    }
# put the reverse and forward reads to own folders
fnFs <- sort(list.files("input_folder", pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files("input_folder", pattern="_R2_001.fastq", full.names = TRUE))

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
#print(sample.names)

# make a variable for the reverse and forward reads for the filterandtrim function
filtFs <- file.path("output_folder", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path("output_folder", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names


out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(truncf,truncr),
              maxN=maxns, maxEE=c("inf","inf"), truncQ=0, rm.phix=TRUE,
              compress=TRUE, multithread=TRUE, verbose=TRUE)
#print(filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(0,0),
           #   maxN=10000000, maxEE=c("Inf","Inf"), truncQ=0, rm.phix=TRUE,
           #   compress=TRUE, multithread=TRUE, verbose=TRUE))
sink("summary.txt")
	cat("\n\n\n")
  cat("Sample name")
	print(out)
	cat("\n\n\n")
sink()
print(list.files("output_folder"))
#file.create("summary.txt")
#write(out,file="summary.txt",append=TRUE)
write.table(out, file ="summary.tsv", row.names = TRUE)

#make a output tar package named contigs.tar and qzip
system("gzip output_folder/*.fq")
system("cd output_folder && tar cf ../filtered.tar *")

#EOF