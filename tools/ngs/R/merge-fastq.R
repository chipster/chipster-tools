# TOOL merge-fastq.R: "Merge FASTQ" (Merges FASTQ files. Files can be gzip compressed.)
# INPUT reads{...}.data: "Fastq files" TYPE GENERIC 
# OUTPUT OPTIONAL merged.fq
# OUTPUT OPTIONAL merged.fq.gz
# PARAMETER order: "Order by file name" TYPE [yes, no] DEFAULT yes (Concatenate files in alphabetical order by file name.)
# PARAMETER compress: "Compress merged file with gzip" TYPE [yes, no] DEFAULT yes (Compress merged file with gzip.)

# AMS 2017-03-09

source(file.path(chipster.common.path, "zip-utils.R"))
source(file.path(chipster.common.path, "tool-utils.R"))

# Read input names
input.names <- read.table("chipster-inputs.tsv", header=F, sep="\t")
# Make input list
input.list <- vector(mode="character", length=0)
if (order == "yes"){
	# Note: files are sorted by display names but internal names are put into input.list
	input.names.sorted <- input.names[order(input.names[,2]),]
	for (i in 1:nrow(input.names)) {
		input.list <- c(input.list, paste(input.names.sorted[i,1]))
	}
}else{
	for (i in 1:nrow(input.names)) {
		input.list <- c(input.list, paste(input.names[i,1]))
	}
}

# Go throug input list and concatenate to output.
for (i in 1:length(input.list)) {
	if (isGZipFile(input.list[i])) {
		system(paste("zcat", input.list[i], ">> merged.fq"))
	}else{
		system(paste("cat", input.list[i], ">> merged.fq"))
	}	
}

# Compress output if required.
if (compress == "yes"){
	system(paste("gzip merged.fq"))
}
