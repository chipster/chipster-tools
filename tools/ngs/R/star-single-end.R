# TOOL star-single-end.R: "STAR for single end reads" (STAR single end)
# INPUT reads{...}.fq: "Reads" TYPE GENERIC
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

# Input fastq names
reads1 <- paste(grep("reads", input.names[,1], value = TRUE), sep="", collapse=",")


# command
command <- paste(star.binary, "--genomeDir", path.star.index, "--readFilesIn", reads1, "--outSAMtype BAM SortedByCoordinate", "--twopassMode Basic", "--runThreadN", chipster.threads.max)
# Use GTF if provided
if (fileOk("annotation.gtf")){
	command <- paste(command, "--sjdbGTFfile annotation.gtf")
}

# Run STAR
#system(command)
runExternal(command)

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
