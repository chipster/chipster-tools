# TOOL fastq-filter.R: "Filter sequences based on the number of expected errors" (Discards sequences which have more than the user specified number of expected errors. This tool can be used for example to filter out bad quality contigs, if paired reads were merged with VSEARCH which preserves the base quality information. )
# INPUT contigs.tar: "Tar package containing the FASTQ files" TYPE GENERIC
# OUTPUT filtered_contigs.tar
# OUTPUT summary.tsv
# PARAMETER maxee: "Discard sequences with more than the specified number of expected errors" TYPE DECIMAL FROM 0 DEFAULT 1 (Decimal number to discard bad quality sequences. Decimal from 0 -, Default 1.)
# PARAMETER qmax: "Maximum quality score accepted when reading FASTQ files" TYPE INTEGER FROM 0 DEFAULT 50 (Specify the maximum quality score accepted when reading FASTQ files. For recent Sanger/Illumina 1.8+ files value 41 is usual. For Nanopore data it can be higher.)

# ES 21.7.2021
# ES  11.8.2022

source(file.path(chipster.common.lib.path, "tool-utils.R"))
source(file.path(chipster.common.lib.path, "zip-utils.R"))

# check out if the file is compressed and if so unzip it
unzipIfGZipFile("reads.tar")

# binary
binary <- c(file.path(chipster.tools.path, "vsearch", "vsearch"))
version <- system(paste(binary, "--version"), intern = TRUE)
documentVersion("Vsearch", version)

# Read the contents of the tar file into a list
system("tar tf contigs.tar > tar.contents")
file.list <- scan("tar.contents", what = "", sep = "\n")

# Check that the input is a valid tar file
if (length(file.list) == 0) {
    stop(paste("CHIPSTER-NOTE: ", "It seems your input file is not a valid Tar package. Please check your input file."))
}
# Check if tar packa contains folders
if (grepl("/", file.list[1])) {
    stop(paste("CHIPSTER-NOTE: ", "It seems your Tar package contains folders. The FASTQ files need to be in the root of the package, not in subfolders."))
}
# Input is a tar file: make input and output folder and create file samples.fastq.txt
system("mkdir input_folder")
system("mkdir output_folder")

# make 4 vectors for summary.tsv dataframe
vector_sequences <- c()
vector_discard <- c()
vector_proportion <- c()


# untar the tar package to input_folder and list the filenames
untar("contigs.tar", exdir = "input_folder")
filenames <- list.files("input_folder", full.names = TRUE)
txt_filenames <- list.files("input_folder")
# file names
sample.names <- sapply(strsplit(basename(txt_filenames), "_"), `[`, 1)

x <- 1
# run the fastq_filter tool for every file
for (file in filenames) {
    # name the input and output file
    output_fastq <- paste0("output_folder/", sample.names[x])
    x <- x + 1
    # make the fastq_filter command
    command <- paste(binary, "--fastq_filter", file, "--fastq_maxee", maxee, "--fastq_qmax", qmax, "--fastqout", output_fastq, ">>summary_test.txt 2>&1")
    # run command
    runExternal(command)
    documentCommand(command)
    system(command)
}
# create summary.tsv dataframe and add filter info for each sample
summary_data <- readLines("summary_test.txt")
all <- 0
disc <- 0
for (row in summary_data) {
    if (grepl("sequences kept", row)) {
        info <- strsplit(row, " ", fixed = TRUE)
        vector_sequences <- c(vector_sequences, info[[1]][1])
        vector_discard <- c(vector_discard, info[[1]][8])
        sum <- strtoi(info[[1]][8]) / (strtoi(info[[1]][8]) + strtoi(info[[1]][1])) * 100
        sum <- round(sum, digits = 4)
        vector_proportion <- c(vector_proportion, sum)
        all <- all + strtoi(info[[1]][1])
        disc <- disc + strtoi(info[[1]][8])
    }
}
sample.names <- append(sample.names, "All", 0)
vector_sequences <- append(vector_sequences, all, 0)
vector_discard <- append(vector_discard, disc, 0)
vector_proportion <- append(vector_proportion, "", 0)
summ.data <- data.frame("Sample" = sample.names, "Sequences_kept" = vector_sequences, "Sequences_discarded" = vector_discard, "Percentage_of_discarded" = vector_proportion)
write.table(summ.data, file = "summary.tsv", row.names = FALSE)
# make a output tar package named contigs.tar and qzip all the files
system("gzip output_folder/*.fq")
system("cd output_folder && tar cf ../filtered_contigs.tar *")
