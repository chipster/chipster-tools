# TOOL hisat2-single-end-with-index-building.R: "HISAT2 for single end reads and own genome" (Aligns single end RNA-seq reads to a user-supplied genome. You need to supply the genome in FASTA format.)
# INPUT reads{...}.fq: "Reads to align" TYPE GENERIC
# INPUT OPTIONAL genome.txt: "Genome to align against" TYPE GENERIC
# OUTPUT OPTIONAL hisat.bam
# OUTPUT OPTIONAL hisat.bam.bai
# OUTPUT OPTIONAL hisat.log
# PARAMETER library.type: "Library type" TYPE [fr-unstranded: fr-unstranded, fr-firststrand: fr-firststrand, fr-secondstrand: fr-secondstrand] DEFAULT fr-unstranded (Which library type to use. For directional\/strand specific library prepartion methods, choose fr-firststrand or fr-secondstrand depending on the preparation method: if the first read \(read1\) maps to the opposite, non-coding strand, choose fr-firststrand. If the first read maps to the coding strand, choose fr-secondstrand. For example for Illumina TruSeq Stranded sample prep, choose fr-firstsrand.)
# PARAMETER OPTIONAL max.multihits: "How many hits to report per read" TYPE INTEGER FROM 1 TO 1000000 DEFAULT 5 (Instructs HISAT2 to report up to this many alignments to the reference for a given read.)
# PARAMETER OPTIONAL quality.format: "Base quality encoding used" TYPE [sanger: "Sanger - Phred+33", phred64: "Phred+64"] DEFAULT sanger (Quality encoding used in the fastq file.)
# PARAMETER OPTIONAL min.intron.length: "Minimum intron length" TYPE INTEGER FROM 1 TO 1000 DEFAULT 20 (Sets minimum intron length. Default: 20)
# PARAMETER OPTIONAL max.intron.length: "Maximum intron length" TYPE INTEGER FROM 1 TO 1000000 DEFAULT 500000 (Sets maximum intron length. Default: 500000)
# PARAMETER OPTIONAL no.softclip: "Disallow soft clipping" TYPE [nosoft: "No soft-clipping", yessoft: "Use soft-clipping"] DEFAULT yessoft (Is soft-cliping used. By default HISAT2 may soft-clip reads near their 5' and 3' ends.)
# PARAMETER OPTIONAL dta: "Require long anchor lengths for subsequent assembly" TYPE [nodta: "Don't require", yesdta: "Require"] DEFAULT nodta (With this option, HISAT2 requires longer anchor lengths for de novo discovery of splice sites. This leads to fewer alignments with short-anchors, which helps transcript assemblers improve significantly in computation and memory usage.)

# AO 30.5.2017 First version
# EK 18.10.2017 Polishing

## Source required functions
source(file.path(chipster.common.path, "zip-utils.R"))
source(file.path(chipster.common.path, "tool-utils.R"))
source(file.path(chipster.common.path, "bam-utils.R"))

## Helper functions
#Unzips a list of files
unzipInputs <- function(names) {
	for (i in 1:nrow(names)) {
		unzipIfGZipFile(names[i,1])	
	}
}

# Echoes command in log file if debug == TRUE
debugPrint <- function(command) {
	if (debug) {
		system(paste("echo ",command, ">> debug.log"))
	}
}

## Options
# Prefer fixed representation over exponential
options(scipen = 10)
# Debug mode, change debug to TRUE or FALSE, depending do you want debug prints or not
debug <- TRUE
debugPrint("")
debugPrint("DEBUG MODE IS ON")

# setting up HISAT binaries (and paths)
hisat.binary <- file.path(chipster.tools.path, "hisat2", "hisat2")
samtools.binary <- file.path(chipster.tools.path, "samtools", "samtools")
hisat.index.binary <- file.path(chipster.tools.path, "hisat2", "hisat2-build")

# Get input name
input.names <- read.table("chipster-inputs.tsv", header=F, sep="\t")

# check out if the file is compressed and if so unzip it
unzipInputs(input.names)

# Do indexing
system("echo Indexing the genome... >> hisat.log")
system("echo >> hisat.log")
index.command<- paste("bash -c '", hisat.index.binary, "-p", chipster.threads.max, "genome.txt", "own_genome", "2>>hisat.log", "'")
system(paste("echo ", index.command, " >> hisat.log"))
system(index.command)


## Parse parameters and store them into hisat.parameters
hisat.parameters <- ""
# Reads to align
# Parse the read names from input files
reads.parsed <- paste(grep("reads", input.names[,1], value = TRUE), sep="", collapse=",")
hisat.parameters <-paste(hisat.parameters, "-U", reads.parsed)
# Quality score format
if (quality.format=="phred64") {
	hisat.parameters <- paste(hisat.parameters, "--phred64")
} else {
	hisat.parameters <- paste(hisat.parameters, "--phred33") 
}
# Intron length, defaluts are 20 and 500 000
hisat.parameters <- paste(hisat.parameters, "--min-intronlen", min.intron.length)
hisat.parameters <- paste(hisat.parameters, "--max-intronlen", max.intron.length)
# Specify strand-specific information: the default is unstranded
if (library.type == "fr-firststrand") {
	hisat.parameters <- paste(hisat.parameters, "--rna-strandness F")
} else if (library.type == "fr-secondstrand") {
	hisat.parameters <- paste(hisat.parameters, "--rna-strandness R")	
}
# Organism
# -x declares the basename of the index for reference genome
hisat.parameters <- paste(hisat.parameters, "-x", "own_genome")
# Set environment variable that defines where indexes locate, HISAT2 requires this
#Sys.setenv(HISAT2_INDEXES = "/opt/chipster/tools/genomes/indexes/hisat2")
# We have created own genome -> no need fot environment variable
# Known splice sites
if (file.exists("splicesites.txt")) {
    hisat.parameters <- paste(hisat.parameters, "--known-splicesite-infile", "splicesites.txt")
}
# How many hits is a read allowed to have
hisat.parameters <-paste(hisat.parameters, "-k", max.multihits)
# Allow soft-clipping, by default soft-clipping is used
if (no.softclip=="nosoft") {
    hisat.parameters <- paste(hisat.parameters, "--no-softclip")
}
# Is longer anchor lengths required
if (dta=="yesdta") {
    hisat.parameters <- paste(hisat.parameters, "--dta")
}

## Set parameters that are not mutable via Chipster
# Threads that hisat uses
hisat.parameters <- paste(hisat.parameters, "-p" , chipster.threads.max)
# Name of the output file
hisat.parameters <- paste(hisat.parameters,  "-S", "hisat.sam")
# Forward errors to hisat.log
hisat.parameters <- paste(hisat.parameters,  "2>> hisat.log")
# Suppress SAM records for reads that failed to align
hisat.parameters <- paste(hisat.parameters, "--no-unal")

#Print the HISAT2_INDEXES into debug
debugPrint("")
debugPrint("HISAT2_INDEXES:")
debugPrint("$HISAT2_INDEXES")

# Print parameters into log
debugPrint("HISAT PARAMETERS")
debugPrint(toString(hisat.parameters))

## Run HISAT

# Note a single ' at the beginning, it allows us to use special characters like >
command <- paste("bash -c '", hisat.binary)

# Add the parameters
command <- paste(command, hisat.parameters)

# Close the command with a ', because there is a opening ' also
command <- paste(command,"'") 
# Print the command to the hisat.log file
debugPrint(command)

# Run command
system(command)

## Run samtools
# Convert SAM file into BAM file and index bam file
# Parameters:
# 	-b Output in the BAM format.
#	-S Ignored for compatibility with previous samtools versions. Previously this option was required if input was in SAM format, but now the correct
#		format is automatically detected by examining the first few characters of input.
# Convert SAM to BAM
debugPrint("")
debugPrint("SAMTOOLS")
samtools.view.command <- paste(samtools.binary, "view -bS hisat.sam > hisat.tmp.bam")
debugPrint(samtools.view.command)
system(samtools.view.command)
# Index bam, this produces a "hisat.sorted.bam" file
samtools.sort.command <- paste(samtools.binary, "sort hisat.tmp.bam hisat.sorted")
debugPrint(samtools.sort.command)
system(samtools.sort.command)

# Do not return empty BAM files
if (fileOk("hisat.sorted.bam", minsize=100)){
	# Rename result files
	system("mv hisat.sorted.bam hisat.bam")
	# Change file names in BAM header to display names
	displayNamesToBAM("hisat.bam")
	# Index BAM
	system(paste(samtools.binary, "index hisat.bam > hisat.bam.bai"))
}


Sys.unsetenv("HISAT2_INDEX")

# Append the debug.log into hisat.log
system("cat debug.log >> hisat.log")

# Substitute display names to log for clarity
displayNamesToFile("hisat.log")

# Handle output names
#
# read input names
inputnames <- read_input_definitions()

# Determine base name
basename <- strip_name(inputnames$reads001.fq)

# Make a matrix of output names
outputnames <- matrix(NA, nrow=2, ncol=2)
outputnames[1,] <- c("hisat.bam", paste(basename, ".bam", sep =""))
outputnames[2,] <- c("hisat.bam.bai", paste(basename, ".bam.bai", sep =""))

# Write output definitions file
write_output_definitions(outputnames)

#EOF

