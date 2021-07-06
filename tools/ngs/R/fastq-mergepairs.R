# TOOL fastq-mergepairs.R: "Combine paired reads to contigs with Vsearch" (Combines paired reads to sequence contigs within each sample, and puts all the resulting FASTQ files in one tar-package. Input file is a single Tar package containing all the FASTQ files, which can be gzipped. You can make a Tar package of your FASTQ files using the Utilities tool Make a tar package. The tool tries to assign the FASTQ files into samples based on the file names, but you can also provide a file containing this information, please see the manual. This tool is based on the Vsearch fastq_mergepairs command.)
# INPUT reads.tar: "Tar package containing the FASTQ files" TYPE GENERIC
# OUTPUT contigs.tar
# OUTPUT summary.txt
# OUTPUT summary_stats.tsv
# OUTPUT samples.fastqs.txt
# PARAMETER maxdiff: "Maximun difference" TYPE INTEGER FROM 0 TO 100 DEFAULT 10 (The Maximum number of non-matching nucleotides allowed in the overlap region 0 - 100, default 10)
# PARAMETER maxdiffpct: "Maximum percentage of non-matching nucleotides" TYPE INTEGER FROM 0 TO 100 DEFAULT 100 (the maximum percentage of non-matching nucleotides allowed in the overlap region. The default value is 100.0%)


# ES 05.07.2021

source(file.path(chipster.common.path,"tool-utils.R"))
source(file.path(chipster.common.path,"zip-utils.R"))

#check out if the file is compressed and if so unzip it
unzipIfGZipFile("reads.tar")

# binary
binary <- c(file.path(chipster.tools.path,"vsearch","vsearch"))
version <- system(paste(binary,"--version"),intern = TRUE)
documentVersion("Vsearch",version)

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


# Input is a tar file: make input and output folder and create file samples.fastq.txt
system("mkdir input_folder")
system("mkdir output_folder")
file.create("samples.fastqs.txt")

# untar the tar package to input_folder and list the filenames
untar("reads.tar", exdir = "input_folder")
filenames <- list.files("input_folder")

# check if the lenght of files in the input_folder is even (all the fastq files have a pair), else error
number <- length(filenames)%%2
if (number != 0 && length(filenames)<2){
    stop(paste('CHIPSTER-NOTE: ',"It seems that some of your fastq files doesn`t have a pair"))
    }

# Sort the filenames _R1 _R2, samples have different names
filenames <- sort(filenames)

#make 2 empty vectors for summary.tsv
vector_names <- c()
vector_info <- c()

x<-1
# run fastq_mergepair for each sample, take _R1 and _R2 files 
while (x < length(filenames)){
    # take the sequence names out
    name <- filenames[x]
    names <- strsplit(name, "_")[[1]]
    # make a output- and input fastq file
    output_fastq <- paste("output_folder/", names, ".fq", sep="")
    input_fastq1 <- paste("input_folder/",name, sep="")
    x = x+1
    # read the reverse read from the fastq pair
    name2 = filenames[x]
    input_fastq2 <- paste("input_folder/",name2, sep="")
    line2 <- paste(names[1]," ", name," ",name2)
    
    # command with parameters
    command <- paste(binary,"--fastq_mergepairs",input_fastq1,"--reverse",input_fastq2,"--eeout",
        "--fastq_maxdiffs",maxdiff,"--fastq_maxdiffpct",maxdiffpct,"--fastq_maxdiffpct 100",
        "--fastqout", output_fastq ,"--label_suffix",names,">>summary.txt 2>&1")

    #add filename to summary.txt file 
    write("\n",file="summary.txt",append=TRUE)
    line = paste("#",names[1])
    write(line,file="summary.txt",append=TRUE)

    # write samples.fastqs.txt file, to check if right samples have been merged
    write(line2,file="samples.fastqs.txt",append=TRUE)
    write("\n",file="summary.txt",append=TRUE)

    #take the sample name to vector
    vector_names <- c(vector_names, names[1])

    # run command
    runExternal(command)
    documentCommand(command)
    system(command)
    x= x+1
    
}
# collect the not merged % from the summary.txt and put it to vector
summary_data <- readLines("summary.txt")
for (row in summary_data){
    if (grepl("Not",row)){
        count <- strsplit(row, " Not merged (",fixed=TRUE)
        count <- count[[1]][2]
        count <- strsplit(count, "%)",fixed=TRUE)
        vector_info <- c(vector_info,count[[1]][1])
    }
}
# make a data.frame and summary.tsv file, where are the sample IDs and not_merged info
summ.data <- data.frame("Sample"= vector_names, "Not_merged"= vector_info)
write.table(summ.data, file ="summary_stats.tsv", row.names = FALSE)

#make a output tar package named contigs.tar and qzip
system("gzip output_folder/*.fq")
system("cd output_folder && tar cf ../contigs.tar *")



