# TOOL picard-merge-bam-alignment.R: "Merge aligned and unaligned BAM" (Merge sorted BAM alignment and unaligned, tagged BAM file. Make sure the input files are assigned correctly!)
# INPUT unmapped.bam: "Unaligned BAM" TYPE BAM
# INPUT aligned.bam: "Aligned BAM" TYPE BAM
# INPUT OPTIONAL own.gtf: "Own GTF reference file" TYPE GENERIC
# OUTPUT merged_tagged.bam
# OUTPUT merged_tagged.bam.bai
# PARAMETER reference: "Reference genome" TYPE ["FILES genomes/fasta .fa"] DEFAULT "SYMLINK_TARGET genomes/fasta/default .fa" (Use same reference as in the alignment!)
# PARAMETER OPTIONAL organism: "Reference GTF" TYPE [other, "FILES genomes/gtf .gtf"] DEFAULT other (GTF file to be used in tagging. No need to select anything here if you are using your own GTF file.)


# 2016-10-31 ML
# 2017-05-08 ML add sort BAM to this tool
# 2017-07-19 ML Naming of outputs according to the input names
# 2018-10-10 ML Combine with Tag reads with gene names tool

# For testing, not run:
# OUTPUT OPTIONAL log.txt
# OUTPUT OPTIONAL stderr.log
# OUTPUT OPTIONAL hs_err_pid17747.log

## Source required functions
source(file.path(chipster.common.path, "tool-utils.R"))
# runExternal writes log to stderr.log 

picard.binary <- file.path(chipster.tools.path, "picard-tools", "picard.jar")
#path.reference <- "/opt/chipster/genomes/fasta"
path.reference <- file.path(chipster.tools.path, "genomes", "fasta")
path.dropseq <- c(file.path(chipster.tools.path, "drop-seq_tools"))

# create symlink
command <- paste("ln -s ", path.reference, "/", reference, ".fa ", reference, ".fa", sep="" )
runExternal(command, checkexit = FALSE)		

# create dictionary (dict) (takes less than minute)
command <- paste("java -Xmx2g -jar ", picard.binary, " CreateSequenceDictionary R=", reference, ".fa O=" ,reference, ".dict", sep="")
runExternal(command, checkexit = FALSE)

# Sort BAM (Picard):
command <- paste("java -Xmx2g -jar", picard.binary, "SortSam I=aligned.bam O=aligned_sorted.bam SO=queryname")
runExternal(command, checkexit = FALSE)

# Merge files (Picard):
command <- paste("java -Xmx2g -jar ", picard.binary, " MergeBamAlignment UNMAPPED_BAM=unmapped.bam ALIGNED_BAM=aligned_sorted.bam O=merged.bam R=" ,reference, ".fa INCLUDE_SECONDARY_ALIGNMENTS=false PAIRED_RUN=false", sep="") 
runExternal(command, checkexit = FALSE)

# Only index if BAM not empty to prevent returning an empty .bai file
if (fileOk("merged.bam", minsize=100)){
	# Index BAM
	samtools.binary <- file.path(chipster.tools.path, "samtools", "samtools")
	system(paste(samtools.binary, "index merged.bam > merged.bam.bai"))
}else{
	system("mv stderr.log merge.log")
}


## Annotate/tag reads with gene names:

# check out if the file is compressed and if so unzip it
source(file.path(chipster.common.path, "zip-utils.R"))
unzipIfGZipFile("own.gtf")

# command start
# If using own GTF:
if (file.exists("own.gtf")){
	library(tools)
	dir <- getwd()
	dir <- file_path_as_absolute(dir)
	command <- paste(path.dropseq, "/TagReadWithGeneExon I=merged.bam O=merged_tagged.bam ANNOTATIONS_FILE=", dir, "/own.gtf TAG=GE", sep="") 
	
# if using one of the GTFs available on Chipster:			
}else{
	#gtf.path <- "/opt/chipster/genomes/gtf/"
	file.path(chipster.tools.path, "genomes", "gtf")
	command <- paste(path.dropseq, "/TagReadWithGeneExon I=merged.bam O=merged_tagged.bam ANNOTATIONS_FILE=", gtf.path, organism, ".gtf TAG=GE", sep="")
}

# run the tool
runExternal(command, checkexit = FALSE)

# Only index if BAM not empty to prevent returning an empty .bai file
if (fileOk("merged_tagged.bam", minsize=100)){
	# Index BAM
	samtools.binary <- file.path(chipster.tools.path, "samtools", "samtools")
	system(paste(samtools.binary, "index merged_tagged.bam > merged_tagged.bam.bai"))
}


# Handle output names
# Source read_input_definitions and strip_name functions
source(file.path(chipster.common.path, "tool-utils.R"))
# read input names and strip file extension
inputnames <- read_input_definitions()
input1name <- inputnames$unmapped.bam # name from the unmapped.bam input
input1namestripped <-strip_name(input1name)
#write the input file name into log
write(input1namestripped, file = "merge.log")
#
## Make a matrix of output names
## These override the default ones
outputnames <- matrix(NA, nrow=2, ncol=2)
outputnames[1,] <- c("merged_tagged.bam", paste(input1namestripped, "_merged_tagged.bam", sep = ""))
outputnames[2,] <- c("merged_tagged.bam.bai", paste(input1namestripped, "_merged_tagged.bam.bai", sep = ""))
#outputnames[3,] <- c("merged_tagged.bam", paste(input1namestripped, "_merged_tagged.bam", sep = ""))
#outputnames[4,] <- c("merged_tagged.bam.bai", paste(input1namestripped, "_merged_tagged.bam.bai", sep = ""))
## Write output definitions file
write_output_definitions(outputnames)

# system("ls -l >> log.txt")

# EOF