# TOOL dada2-merge.R: "Make contigs and remove chimeras" (Giving the filtered fastq files in a tar package and two dada-class objects produced with dada function, this tool combines the files to contigs, makes an asv-table and removes chimeras. For more information check the manual)
# INPUT forward.Rda: "DADA-class object created from forward reads" TYPE GENERIC (dada-class object saved as .Rda file and created with dada function.)
# INPUT reverse.Rda: "DADA-class object created from reverse reads" TYPE GENERIC (dada-class object saved as .Rda file and created with dada function.)
# INPUT reads.tar: "Tar package containing the FASTQ files" TYPE GENERIC
# OUTPUT OPTIONAL seqtab_nochim.Rda
# OUTPUT OPTIONAL log4.txt
# OUTPUT OPTIONAL log.txt
# OUTPUT OPTIONAL preprocessing_summary.tsv
# OUTPUT OPTIONAL sequence_table.tsv
# RUNTIME R-4.1.1

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

# Make input folder
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


# not working!!
sink(file="log4.txt")
  sink.number(type = "output")
    #cat("\nmergePairs:\n")
      #errF <- learnErrors(fnFs, multithread=TRUE)
    mergers <- mergePairs(dadaFs, fnFs, dadaRs, fnRs, verbose=TRUE)
    print(head(mergers[[1]]))
sink()

# Construct sequence table from mergers -> ASV-table Construct a sample-by-sequence observation matrix.
seqtab <- makeSequenceTable(mergers)

#rename the rows by splitting the name with _
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
rownames(seqtab) <- sample.names

# Write a log/summary file
sink(file="log.txt")
    cat("\n")
    cat("After mergePairs command sequence table consist of:\n")
    cat(length(rownames(seqtab)))
    cat(" samples and ")
    cat(length(colnames(seqtab)))
    cat(" amplicon sequence variants\n")
    cat("\nDistribution of sequence lengths:\n")
    table(nchar(getSequences(seqtab)))
    cat("\n")

# run isbimeradenovo / remove chimeras
    seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
    num <- length(colnames(seqtab))-length(colnames(seqtab.nochim))
    cat("Identified ")
    cat(num)
    cat(" bimeras out of ")
    cat(length(colnames(seqtab)))
    cat(" input sequences\n")
    cat("Total amount of ASVs is: ")
    cat(length(colnames(seqtab.nochim)))
    cat("\n")
sink()

# track reads through the pipeline and make a tsv table
getN <- function(x) sum(getUniques(x))
track <- cbind(sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
write.table(track, file="preprocessing_summary.tsv", sep="\t", row.names=TRUE, col.names=T, quote=F)


# print out asv sequence table rename asv:s
# rename sequences to asv1, asv2... makes it easier
seqtab.nochim2 <- seqtab.nochim
names <- c()
x=0
while (x<length(colnames(seqtab.nochim2))){
    new = paste("asv",x,sep="")
    names <- c(names, new)
    x = x+1
}
colnames(seqtab.nochim2) <- names
write.table(seqtab.nochim2, file="sequence_table.tsv", sep="\t", row.names=TRUE, col.names=T, quote=F)

# save the object as .Rda 
save(seqtab.nochim, file = "seqtab_nochim.Rda")