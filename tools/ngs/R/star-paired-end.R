# TOOL star-paired-end.R: "STAR for paired end reads" (STAR)
# INPUT reads{...}.fq: "Reads" TYPE GENERIC
# INPUT OPTIONAL reads1.txt: "List of read 1 files" TYPE GENERIC
# INPUT OPTIONAL reads2.txt: "List of read 2 files" TYPE GENERIC
# INPUT OPTIONAL annotation.gtf: "Optional GTF file" TYPE GENERIC
# OUTPUT OPTIONAL alignment.bam
# OUTPUT OPTIONAL alignment.bam.bai
# OUTPUT OPTIONAL Log_progress.txt
# OUTPUT OPTIONAL Log_final.txt
# PARAMETER organism: "Genome" TYPE [Homo_sapiens.GRCh38_20, Homo_sapiens.GRCh38_20_gtf] DEFAULT Homo_sapiens.GRCh38_20_gtf (Genome or transcriptome that you would like to align your reads against.)

source(file.path(chipster.common.path, "tool-utils.R"))

# check out if the file is compressed and if so unzip it
source(file.path(chipster.common.path, "zip-utils.R"))
input.names <- read.table("chipster-inputs.tsv", header=F, sep="\t")
for (i in 1:nrow(input.names)) {
	unzipIfGZipFile(input.names[i,1])	
}

# setting up STAR
star.binary <- c(file.path(chipster.tools.path, "STAR", "STAR"))
path.star.index <- c(file.path(chipster.tools.path, "genomes", "indexes", "star", organism))
samtools.binary <- c(file.path(chipster.tools.path, "samtools", "samtools"))

# Input files
if (fileOk("reads1.txt", 0) && fileOk("reads2.txt", 0)){
	# Case: list files exist
	reads1.list <- make_input_list("reads1.txt")
	reads2.list <- make_input_list("reads2.txt")
	if (identical(intersect(reads1.list, reads2.list), character(0))){
		reads1 <- paste(reads1.list, sep="", collapse=",")
		reads2 <- paste(reads2.list, sep="", collapse=",")
	}else{
		stop(paste('CHIPSTER-NOTE: ', "One or more files is listed in both lists."))
	}
}else if (fileOk("reads002.fq") && fileNotOk("reads003.fq")){
	# Case: no list file, but only two fastq inputs
	in.sorted <- input.names[order(input.names[,2]),]
	reads <- grep("reads", in.sorted[,1], value = TRUE)
	reads1 <- reads[1]
	reads2 <- reads[2]
}else{
	# Case: no list files, more than two fastq inputs
	stop(paste('CHIPSTER-NOTE: ', "List file is missing. You need to provide a list of read files for both directions."))
}

# command
command <- paste(star.binary, "--genomeDir", path.star.index, "--readFilesIn", reads1, reads2, "--outSAMtype BAM SortedByCoordinate", "--twopassMode Basic", "--runThreadN", chipster.threads.max)
# Use GTF if provided
if (fileOk("annotation.gtf")){
	command <- paste(command, "--sjdbGTFfile annotation.gtf")
}

# Run STAR
system(command)

# rename result files
system("mv Log.progress.out Log_progress.txt")
system("mv Log.final.out Log_final.txt")
system("mv Aligned.sortedByCoord.out.bam alignment.bam")

# index bam
system(paste(samtools.binary, "index alignment.bam"))

# Determine base name
inputnames <- read_input_definitions()
basename <- strip_name(inputnames$reads001.fq)

# Make a matrix of output names
outputnames <- matrix(NA, nrow=2, ncol=2)
outputnames[1,] <- c("alignment.bam", paste(basename, ".bam", sep =""))
outputnames[2,] <- c("alignment.bam.bai", paste(basename, ".bam.bai", sep =""))

# Write output definitions file
write_output_definitions(outputnames)
