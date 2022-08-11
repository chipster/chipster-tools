# TOOL fastq-mergepairs.R: "Combine paired reads to contigs with VSEARCH" (Combines paired reads to sequence contigs within each sample, and puts all the resulting FASTQ files in one Tar package. Input file is a single Tar package containing all the FASTQ files, which can be gzipped. You can make a Tar package of your FASTQ files using the Utilities tool Make a tar package. The tool tries to assign the FASTQ files into samples based on the file names, but you can also provide a file containing this information, please see the manual. This tool is based on the VSEARCH fastq_mergepairs command.)
# INPUT reads.tar: "Tar package containing the FASTQ files" TYPE GENERIC
# INPUT OPTIONAL input_list.txt: "List of FASTQ files by sample" TYPE GENERIC
# OUTPUT contigs.tar
# OUTPUT summary.txt
# OUTPUT summary_stats.tsv
# OUTPUT samples.fastqs.txt
# PARAMETER maxdiff: "Maximun number of non-matching nucleotides" TYPE INTEGER FROM 0 TO 100 DEFAULT 10 (The maximum number of non-matching nucleotides allowed in the overlap region 0 - 100, the default value is 10)
# PARAMETER maxdiffpct: "Maximum percentage of non-matching nucleotides" TYPE INTEGER FROM 0 TO 100 DEFAULT 100 (The maximum percentage of non-matching nucleotides allowed in the overlap region. The default value is 100.0%)
# PARAMETER OPTIONAL maxns: "Discard input sequences with more than this number of Ns." TYPE INTEGER FROM 0 (Discard input sequences which contain more than the specified number of ambigious bases) 

# ES 05.07.2021
# ES 08.07.2022 # made it look much nicer

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

#make 6 empty vectors for summary.tsv
sample.names <- c()
vector_info <- c()
vector_mergeinfo <- c()
vector_pairs <- c()
vector_merged <- c()
vector_notmerged <- c()

#list the files full names with input_folder/
filenames <- list.files("input_folder", full.names=TRUE)
#if input list selected use it, and make a new list of filenames
if (fileOk("input_list.txt")){
    txt_filenames <- c()
    input <- readLines("input_list.txt")
    for (row in input){
        sample <- strsplit(row,'\t',fixed=TRUE)
        sample.names <- c(sample.names, trimws(sample[[1]][1])) 
        first = trimws(sample[[1]][2])
        second = trimws(sample[[1]][3])
        txt_filenames <- c(txt_filenames,first)
        txt_filenames <- c(txt_filenames,second)
        }
    # if the input list has a different amount of files than the tar package
    if (length(txt_filenames) != length(filenames)){
        stop(paste('CHIPSTER-NOTE: ',"It seems that the list of FASTQ files .txt file has different amount of filenames than the .tar package. Please check manual."))
    }else{ #add the full name
        filenames <- paste0("input_folder/",txt_filenames)
    }
}else{
    # Sort the filenames _R1 _R2, samples have different names
    filenames <- sort(filenames)
    txt_filenames <- list.files("input_folder") #names without the folder name
}

# check if the lenght of files in the input_folder is even (all the fastq files have a pair), else error
number <- length(filenames)%%2
if (number != 0){
    stop(paste('CHIPSTER-NOTE: ',"It seems that some of your fastq files doesn`t have a pair"))
    }

# put the forward files to fnFs, assume that forward reads have the same name but different tag than reverse reads
# forward reads should be before reverse reads after sort() function
forward <- seq(1,length(filenames),by=2)
fnFs <- filenames[forward]

reverse <- seq(2,length(filenames), by=2)
fnRs <- filenames[reverse]

# take out the sample names splitting with _, if input
if (!fileOk("input_list.txt")){
  sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

}

x<-1
y<-1
# run fastq_mergepair for each sample, take _R1 and _R2 files 
for (name in sample.names){
    # make a output- and input fastq file
    output_fastq <- paste0("output_folder/", name, ".fq")
    # for  samples.fastqs.txt
    line2 <- paste(name,'\t',txt_filenames[y],'\t',txt_filenames[y+1])
    #if maxns is selected
    if (!is.na(maxns)){ # command with maxns parameter
        command <- paste(binary,"--fastq_mergepairs",fnFs[x],"--reverse",fnRs[x],"--eeout",
            "--fastq_maxdiffs",maxdiff,"--fastq_maxdiffpct",maxdiffpct, "--fastq_maxns",maxns,
            "--fastqout", output_fastq ,"--label_suffix",name,">>summary.txt 2>&1") 
    }else{# else without --fastq_maxns
    # command without maxns parameter
    command <- paste(binary,"--fastq_mergepairs",fnFs[x],"--reverse",fnRs[x],"--eeout",
        "--fastq_maxdiffs",maxdiff,"--fastq_maxdiffpct",maxdiffpct,
        "--fastqout", output_fastq ,"--label_suffix",name,">>summary.txt 2>&1")
    }
    x<-x+1
    y<-y+2

    #add filename to summary.txt file 
    write("\n",file="summary.txt",append=TRUE)
    line = paste("#",name)
    write(line,file="summary.txt",append=TRUE)
    # write samples.fastqs.txt file, to check if right samples have been merged
    write(line2,file="samples.fastqs.txt",append=TRUE)
    write("\n",file="summary.txt",append=TRUE)

    # run command
    runExternal(command)
    documentCommand(command)
    system(command)
    write("--------------------------------------------------------------",file="summary.txt",append=TRUE)
}

# make vectors for summary_stats.tsv file
summary_data <- readLines("summary.txt")
for (row in summary_data){
    #collect pairs
    if (grepl("  Pairs",row)){
        count <- strsplit(row, "  Pairs",fixed=TRUE)
        count<-count[[1]][1]
        count <- strsplit(count," ")
        len <- length(count[[1]])
        vector_pairs<- c(vector_pairs,count[[1]][len])
    }
    # collect merged
    if (grepl("  Merged",row)){
        count <- strsplit(row, "  Merged",fixed=TRUE)
        count<-count[[1]][1]
        count <- strsplit(count," ")
        len <- length(count[[1]])
        vector_merged<- c(vector_merged,count[[1]][len])
    }

    # collect the not merged and % from the summary.txt and put it to vector
    if (grepl("Not",row)){
        count <- strsplit(row, "  Not merged (",fixed=TRUE)
        count2 <- count[[1]][1]
        count2 <- strsplit(count2," ")
        len <- length(count2[[1]])
        vector_notmerged<- c(vector_notmerged,count2[[1]][len])
        count <- count[[1]][2]
        count <- strsplit(count, "%)",fixed=TRUE)
        vector_info <- c(vector_info,count[[1]][1])
    }
    #collect mean expected error
    if (grepl("merged sequences",row)){
        count <- strsplit(row, " Mean",fixed=TRUE)
        count <- count[[1]][1]
        count <- strsplit(count, "%",fixed=TRUE)
        vector_mergeinfo <- c(vector_mergeinfo,count[[1]][1])
    }
}

# make a data.frame and summary.tsv file, where are the sample IDs and not_merged info
summ.data <- data.frame("Sample"= sample.names,"Pairs"=vector_pairs, "Merged"=vector_merged, "Not_merged"=vector_notmerged, "Percentage_not_merged"= vector_info, "Mean_expected_error" = vector_mergeinfo)
write.table(summ.data, file ="summary_stats.tsv", row.names = FALSE)

#make a output tar package named contigs.tar and qzip
system("gzip output_folder/*.fq")
system("cd output_folder && tar cf ../contigs.tar *")



