# TOOL mothur-makegroup.R: "Combine FASTQ or FASTA files and make a group file" (Combines FASTQ or FASTA files of all samples to one file and creates a Mothur group file for it. Input files can be gzipped. This tool is based on the Mothur tools make.group, merge.files and summary.seqs.)
# INPUT reads{...}: "read files" TYPE GENERIC
# OUTPUT sequences.fasta.gz
# OUTPUT sequences.groups
# OUTPUT sequences-summary.tsv
# PARAMETER input.type: "Input type" TYPE [FASTQ, FASTA] DEFAULT FASTQ (Input file type.)


# EK 6.4.2021
# OUTPUT log1.txt


source(file.path(chipster.common.path, "tool-utils.R"))
source(file.path(chipster.common.path, "zip-utils.R"))

# check out if the file is compressed and if so unzip it
input.names <- read.table("chipster-inputs.tsv", header=F, sep="\t")
for (i in 1:nrow(input.names)) {
	unzipIfGZipFile(input.names[i,1])	
}

# Convert to FASTA if necessary
# Binary
binary <- c(file.path(chipster.tools.path, "fastx", "bin", "fastq_to_fasta"))
# Go through the list and convert or rename depending on input type
for (i in 1:nrow(input.names)) {
  newname <- paste(input.names[i,1], ".fasta", sep = "")
  if (input.type == "FASTQ"){
    command <- paste(binary, "-n -i", input.names[i,1],"-o", newname)
    system(command)
  }else{
    file.rename(paste(input.names[i,1]), newname)
  }
}
# Binary
binary <- c(file.path(chipster.tools.path,"mothur","mothur"))
version <- system(paste(binary,"--version"),intern = TRUE)
documentVersion("Mothur",version)

# Collect names of the input files
input.names <- read.table("chipster-inputs.tsv", header=F, sep="\t")

# Set names for first input
fasta_names <- paste(strip_name(paste(input.names[1,1])), ".fasta", sep="")
group_names <- strip_name(paste(input.names[1,2]))

# Loop through the rest of the inputs
for (i in 2:nrow(input.names)) {
  fasta <- paste(strip_name(paste(input.names[i,1])), ".fasta", sep="")
  fasta_names <- paste(fasta_names, fasta, sep="-") 
  group <- strip_name(paste(input.names[i,2]))
  group_names <- paste(group_names, group, sep="-")
}

makegroup.options <- paste("make.group(fasta=",fasta_names, ", groups=",group_names, ")", sep = "")
documentCommand(makegroup.options)
write(makegroup.options,"makegroup.mth",append = FALSE)
command1 <- paste(binary,"makegroup.mth","> log1.txt")
system(command1)

# Rename groups file
filename <- list.files(path =".", pattern="*groups")
system(paste("mv", filename[1], "sequences.groups"))

# Merge fasta files
mergefiles.options <- paste("merge.files(input=", fasta_names, ", output=sequences.fasta)",  sep = "")
documentCommand(mergefiles.options)
write(mergefiles.options,"mergefiles.mth",append = FALSE)
command2 <- paste(binary,"mergefiles.mth","> log2.txt")
system(command2)

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
