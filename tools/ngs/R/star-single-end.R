# TOOL star-single-end.R: "STAR for single end reads and human genome" (Aligns single end RNA-seq reads to a genome. Alignment results are given in a BAM file, which is automatically indexed and hence ready to be viewed in Chipster genome browser. )
# INPUT reads{...}.fq: "Reads" TYPE GENERIC
# INPUT OPTIONAL annotation.gtf: "Optional GTF file" TYPE GENERIC
# OUTPUT OPTIONAL alignment.bam
# OUTPUT OPTIONAL alignment.bam.bai
# OUTPUT OPTIONAL Log_progress.txt
# OUTPUT OPTIONAL Log_final.txt
# PARAMETER organism: "Genome" TYPE [Homo_sapiens.GRCh38] DEFAULT Homo_sapiens.GRCh38 (Genome that you would like to align your reads against.)
# PARAMETER OPTIONAL alignments.per.read: "Maximum alignments per read" TYPE INTEGER DEFAULT 10 (Maximum number of multiple alignments allowed for a read: if exceeded, the read is considered unmapped.)
# PARAMETER OPTIONAL mismatches.per.pair: "Maximum mismatches per read" TYPE INTEGER DEFAULT 10 (Maximum number of mismatches per pair. Use value 999 to switch off this filter.)

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
samtools.binary <- c(file.path(chipster.tools.path, "samtools", "samtools"))

# Input fastq names
reads1 <- paste(grep("reads", input.names[,1], value = TRUE), sep="", collapse=",")


# command
command <- paste(star.binary, "--genomeDir", path.star.index, "--readFilesIn", reads1, "--outSAMtype BAM SortedByCoordinate", "--twopassMode Basic", "--runThreadN", chipster.threads.max, "--alignSJoverhangMin 8", "--alignSJDBoverhangMin 1", "--outFilterType BySJout", "--outFilterMultimapNmax", alignments.per.read, "--outFilterMismatchNmax", mismatches.per.pair)
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
#system(command)
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
