# TOOL dada2-merge.R: "Combine paired reads to contigs with DADA2" (Given the filtered FASTQ files in a tar package and two dada-class objects produced with the tool Sample inference, this tool merges the files to contigs. The tar package needs to be the same as the tar package given to the Sample Inference tool.)
# INPUT forward.Rda: "DADA-class object of forward reads" TYPE GENERIC (dada-class object saved as .Rda file and created with dada function.)
# INPUT reverse.Rda: "DADA-class object of reverse reads" TYPE GENERIC (dada-class object saved as .Rda file and created with dada function.)
# INPUT reads.tar: "Tar package containing the FASTQ files" TYPE GENERIC (Tar package containing those FASTQ files which were used to create the dada-class objects.)
# OUTPUT contigs.Rda
# OUTPUT OPTIONAL contigs_summary.tsv
# PARAMETER minoverlap: "The minimum length of the overlap required for merging the forward and reverse reads" TYPE INTEGER FROM 0 DEFAULT 12 (By default the overlap area should be at least 12 base pairs long.)
# PARAMETER maxmismatch: "The maximum number of mismatches allowed in the overlap region" TYPE INTEGER FROM 0 DEFAULT 0 (By default no mismatches are allowed in the overlap region.)
# RUNTIME R-4.1.1-asv

# ES 11.08.2022
# OUTPUT OPTIONAL contigs.txt
# OUTPUT OPTIONAL contigs2.txt
# PARAMETER OPTIONAL mock: "Name of the mock community if you have co-sequenced a mock community" TYPE STRING (If you have co-sequenced a mock community, you can remove it from the phyloseq object by giving the name of the community as a parameter. Most likely the name is Mock)

source(file.path(chipster.common.path,"tool-utils.R"))
source(file.path(chipster.common.path,"zip-utils.R"))

# load library dada2
library(dada2)
#packageVersion("dada2")

# Load input data 
load("forward.Rda")
load("reverse.Rda")

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

# Make an input folder
system("mkdir input_folder")

# untar the tar package to input_folder and list the filenames
untar("reads.tar", exdir = "input_folder")
filenames <- list.files("input_folder", full.names=TRUE)

 # Sort the filenames _R1 _R2, samples have different names
filenames <- sort(filenames)

# check if the lenght of files in the input_folder is even (all the fastq files have a pair), else error
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

#capture.output(mergers<-mergePairs(dadaFs, fnFs, dadaRs, fnRs, minOverlap=minoverlap, maxMismatch=maxmismatch, verbose=TRUE), file="contigs2.txt")
# sink not working!!
#sink(file="contigs.txt", type=c("output","message"))
 # sink.number(type = "message")
    #cat("\nmergePairs:\n")
   # errF <- learnErrors(fnFs, multithread=TRUE, verbose=TRUE)

mergers <- mergePairs(dadaFs, fnFs, dadaRs, fnRs, minOverlap=minoverlap, maxMismatch=maxmismatch, verbose=TRUE)
#sink()

print(mergers$abundance)


save(mergers, file="contigs.Rda" )

#rename the rows by splitting the name with _
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)


# track reads through these tools and make a tsv table
getN <- function(x) sum(getUniques(x))
track <- cbind(sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN))
colnames(track) <- c("Forward dada object", "Reverse dada object", "After make contigs")
rownames(track) <- sample.names
write.table(track, file="contigs_summary.tsv", sep="\t", row.names=TRUE, col.names=T, quote=F)


