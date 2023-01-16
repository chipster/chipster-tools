# TOOL mothur-makecount.R: "Combine FASTQ files into one FASTA file and make a Mothur count file" (Combines FASTQ files of all samples to one FASTA file and creates a Mothur count file for it. Input file should be a tar file and the files in the tar package can be gzipped.)
# INPUT reads.tar: "Tar package containing the FASTQ files" TYPE GENERIC
# OUTPUT sequences.count_table
# OUTPUT sequences-summary.tsv
# OUTPUT sequences.fasta.gz
# RUNTIME R-4.1.1

# ES 1.12.2022 new mothur version needs count 
#Fastx not in new ubuntu, emboss is. Other possibility is mothur fastq.info
# OUTPUT log1.txt
# OUTPUT log2.txt
# OUTPUT log3.txt

source(file.path(chipster.common.path, "tool-utils.R"))
source(file.path(chipster.common.path, "zip-utils.R"))

# Binary emboss
emboss.binary <- c(file.path(chipster.tools.path, "emboss-20.04", "bin", "seqret"))
version <- system(paste(emboss.binary,"--version"),intern = TRUE)

# Go through the input files. Uncompress if necessary. Untar if necessary. 
fasta_names <- ""
group_names <- ""

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

#list the full file names
#filepaths <- list.files("input_folder", full.names=TRUE)
filenames <- list.files("input_folder")

# Go through the files and gzip. Files that are already gzipped will be skipped
for (file in filenames) {
 
    file_path = paste0("input_folder/",file)
    unzipIfGZipFile(file_path)

    basename <- paste(strip_name(paste(file)))
   
    fastaname <- paste("output_folder/",basename, ".fasta", sep = "")
  
    command <- paste(emboss.binary, "-sequence",file_path,"-outseq", fastaname)
    runExternal(command)   
    fasta_names <- paste(fasta_names, fastaname, sep="-")
    group_names <- paste(group_names, basename, sep="-") 
    }    


# List files now have extra "-" at beginning that has to be removed
fasta_names <- sub(".","", fasta_names)
group_names <- sub(".","",group_names)

#print(fasta_names)
#print(group_names)

# Binary
binary <- c(file.path(chipster.tools.path,"mothur","mothur"))
#binary <- c(file.path(chipster.tools.path,"mothur-1.44.3","mothur"))
version <- system(paste(binary,"--version"),intern = TRUE)
documentVersion("Mothur",version)
makegroup.options <- paste("make.count(fasta=",fasta_names, ", groups=",group_names, ")", sep = "")
documentCommand(makegroup.options)
write(makegroup.options,"makegroup.mth",append = FALSE)
command1 <- paste(binary,"makegroup.mth","> log1.txt")
system(command1)

countseqs.options <- paste("count.seqs(count=output_folder/merge.count_table, compress=f)",sep="") 
documentCommand(countseqs.options)
write(countseqs.options,"batch.mth",append = FALSE)
command <- paste(binary,"batch.mth",">> log3.txt")
system(command)
# rename the result file
system("mv output_folder/merge.full.count_table sequences.count_table")

# Merge fasta files
mergefiles.options <- paste("merge.files(input=", fasta_names, ", output=sequences.fasta)",  sep = "")
documentCommand(mergefiles.options)
write(mergefiles.options,"mergefiles.mth",append = FALSE)
command2 <- paste(binary,"mergefiles.mth","> log2.txt")
system(command2)

system("mv output_folder/sequences.fasta sequences.fasta")

# Create summary file
summaryseqs.options <- paste("summary.seqs(fasta=sequences.fasta")
summaryseqs.options <- paste(summaryseqs.options,", processors=",chipster.threads.max,")",sep = "")
documentCommand(summaryseqs.options)
write(summaryseqs.options,"summary.mth",append = FALSE)
command <- paste(binary,"summary.mth","> log_raw.txt")
system(command)
# Post process output
system("grep -A 10 Start log_raw.txt > summary.tsv")
# Remove one tab to get the column naming look nice:
system("sed 's/^		/	/' summary.tsv > sequences-summary.tsv")

# Gzip output fasta
system("gzip sequences.fasta")

#Sys.sleep(120)