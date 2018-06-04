# TOOL star-paired-end.R: "STAR for paired end reads and human genome" (Aligns paired end RNA-seq reads to a genome. If you have just one pair of read files, Chipster sets reads 1 file and reads 2 file based on file names. If you have more pairs of read files for one sample, you need to provide a list of filenames of the FASTQ files for each direction \(e.g. 1files.txt and 2files.txt\). You can generate the lists with the tool \"Utilities \\\ Make a list of filenames\". Alignment results are given in a BAM file, which is automatically indexed and hence ready to be viewed in Chipster genome browser.)
# INPUT reads{...}.fq: "Reads" TYPE GENERIC
# INPUT OPTIONAL reads1.txt: "List of read 1 files" TYPE GENERIC
# INPUT OPTIONAL reads2.txt: "List of read 2 files" TYPE GENERIC
# INPUT OPTIONAL annotation.gtf: "Optional GTF file" TYPE GENERIC
# OUTPUT OPTIONAL alignment.bam
# OUTPUT OPTIONAL alignment.bam.bai
# OUTPUT OPTIONAL Log_progress.txt
# OUTPUT OPTIONAL Log_final.txt
# PARAMETER organism: "Genome" TYPE [Homo_sapiens.GRCh38.92] DEFAULT Homo_sapiens.GRCh38.92 (Genome that you would like to align your reads against.)
# PARAMETER OPTIONAL alignments.per.read: "Maximum alignments per read" TYPE INTEGER DEFAULT 10 (Maximum number of multiple alignments allowed for a read: if exceeded, the read is considered unmapped.)
# PARAMETER OPTIONAL mismatches.per.pair: "Maximum mismatches per read pair" TYPE INTEGER DEFAULT 10 (Maximum number of mismatches per pair. Use value 999 to switch off this filter.)

source(file.path(chipster.common.path, "tool-utils.R"))
source(file.path(chipster.common.path, "bam-utils.R"))
source(file.path(chipster.common.path, "zip-utils.R"))

# check out if the file is compressed and if so unzip it
input.names <- read.table("chipster-inputs.tsv", header=F, sep="\t")
for (i in 1:nrow(input.names)) {
	unzipIfGZipFile(input.names[i,1])	
}

# setting up STAR
star.binary <- c(file.path(chipster.tools.path, "STAR", "STAR"))
path.star.index <- c(file.path(chipster.tools.path, "genomes", "indexes", "star", organism))
path.gtf <- c(file.path(chipster.tools.path, "genomes", "gtf", organism))
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
command <- paste(star.binary, "--genomeDir", path.star.index, "--readFilesIn", reads1, reads2, "--outSAMtype BAM SortedByCoordinate", "--twopassMode Basic", "--runThreadN", chipster.threads.max, "--alignSJoverhangMin 8", "--alignSJDBoverhangMin 1", "--outFilterType BySJout", "--outFilterMultimapNmax", alignments.per.read, "--outFilterMismatchNmax", mismatches.per.pair)
# Use GTF if provided
if (fileOk("annotation.gtf")){
	command <- paste(command, "--sjdbGTFfile annotation.gtf")
}else{
	gtf.path <- c(file.path(chipster.tools.path, "genomes", "gtf"))
	file.list <- list.files(gtf.path, pattern="\\.gtf$")
	gtf.file <- grep(organism, file.list, value = TRUE)
	gtf.path <- c(file.path(gtf.path, gtf.file[1]))
	command <- paste(command, "--sjdbGTFfile", gtf.path)
}

# Run STAR
runExternal(command)

# rename result files
system("mv Log.progress.out Log_progress.txt")
system("mv Log.final.out Log_final.txt")
system("mv Aligned.sortedByCoord.out.bam alignment.bam")

# Change file named in BAM header to display names
displayNamesToBAM("alignment.bam")

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
