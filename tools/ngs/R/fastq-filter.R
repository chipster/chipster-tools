# TOOL fastq-filter.R: "Filter sequences based on the number of expected errors" (Discards sequences which have more than the user specified number of expected errors. This tool can be used for example to filter out bad quality contigs, if paired reads were merged with VSEARCH which preserves the base quality information. )
# INPUT contigs.tar: "Tar package containing the contig FASTQ files" TYPE GENERIC
# OUTPUT filtered_contigs.tar
# OUTPUT summary.tsv
# PARAMETER maxee: "Discard sequences with more than the specified number of expected errors" TYPE DECIMAL FROM 0 DEFAULT 1 (Decimal number to discard bad quality sequences. Decimal from 0 -, Default 1.)

# ES 21.7.2021

source(file.path(chipster.common.path,"tool-utils.R"))
source(file.path(chipster.common.path,"zip-utils.R"))

#check out if the file is compressed and if so unzip it
unzipIfGZipFile("reads.tar")

# binary
binary <- c(file.path(chipster.tools.path,"vsearch","vsearch"))
version <- system(paste(binary,"--version"),intern = TRUE)
documentVersion("Vsearch",version)

# Read the contents of the tar file into a list
system("tar tf contigs.tar > tar.contents")
file.list <- scan("tar.contents",what = "",sep = "\n")




# Check that the input is a valid tar file
if (length(file.list) == 0) {
    stop(paste('CHIPSTER-NOTE: ',"It seems your input file is not a valid Tar package. Please check your input file."))
}
# Check if tar packa contains folders
if (grepl("/",file.list[1])) {
    stop(paste('CHIPSTER-NOTE: ',"It seems your Tar package contains folders. The FASTQ files need to be in the root of the package, not in subfolders."))
}
# Input is a tar file: make input and output folder and create file samples.fastq.txt
system("mkdir input_folder")
system("mkdir output_folder")

# make 4 vectors for summary.tsv dataframe
vector_names <- c()
vector_sequences <- c()
vector_discard <- c()
vector_proportion <- c()



# untar the tar package to input_folder and list the filenames
untar("contigs.tar", exdir = "input_folder")
filenames <- list.files("input_folder")

# run the fastq_filter tool for every file
for (file in filenames) {
    #take the file name out without .gz
    name <- strsplit(file, '.',fixed=TRUE)
    name <- name[[1]][1]
    vector_names <- c(vector_names,name)
    # name the input and output file
    output_fastq <- paste("output_folder/", name, ".fq", sep="")
    input_fastq <- paste("input_folder/",file, sep="")

    # make the fastq_filter command
    command <- paste(binary,"--fastq_filter", input_fastq, "--fastq_maxee", maxee, "--fastqout", output_fastq, ">>summary_test.txt 2>&1")
    # run command
    runExternal(command)
    documentCommand(command)
    system(command)
}
#create summary.tsv dataframe and add filter info for each sample
summary_data <- readLines("summary_test.txt")
all <- 0
disc <- 0
for (row in summary_data){
    if (grepl("sequences kept",row)){
        info <- strsplit(row, " ",fixed=TRUE)
        vector_sequences <- c(vector_sequences,info[[1]][1])
        vector_discard <- c(vector_discard, info[[1]][8])
        sum <- strtoi(info[[1]][8]) / (strtoi(info[[1]][8])+strtoi(info[[1]][1])) * 100
        sum <- round(sum, digits=4)
        vector_proportion <- c(vector_proportion,sum)
        all <- all + strtoi(info[[1]][1])
        disc <- disc +strtoi(info[[1]][8])
    }
}
vector_names <- append(vector_names,"All",0)
vector_sequences <- append(vector_sequences, all,0)
vector_discard <- append(vector_discard, disc,0)
vector_proportion <- append(vector_proportion,"",0)
summ.data <- data.frame("Sample"= vector_names,"Sequences_kept"=vector_sequences, "Sequences_discarded"= vector_discard, "Percentage_of_discarded" = vector_proportion)
write.table(summ.data, file ="summary.tsv", row.names = FALSE)
#make a output tar package named contigs.tar and qzip all the files
system("gzip output_folder/*.fq")
system("cd output_folder && tar cf ../filtered_contigs.tar *")