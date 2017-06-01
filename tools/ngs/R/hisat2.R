# TOOL hisat2.R: "HISAT2 for single end reads" ( HISAT Aligns single end RNA-seq reads to a genome. Takes only one reads filein ziped format or not )
# INPUT read.fq: "Reads" TYPE GENERIC
# OUTPUT OPTIONAL hisat.bam
# OUTPUT OPTIONAL hisat.log
# PARAMETER organism: "Genome" TYPE ["FILES genomes/indexes/hisat2 .fa"] DEFAULT "SYMLINK_TARGET genomes/indexes/hisat2/default .fa" (Genome or transcriptome that you would like to align your reads against.)
# PARAMETER OPTIONAL quality.format: "Base quality encoding used" TYPE [sanger: "Sanger - Phred+33", phred64: "Phred+64"] DEFAULT sanger (Quality encoding used in the fastq file.)
# PARAMETER OPTIONAL min.intron.length: "Minimum intron length" TYPE INTEGER FROM 1 TO 1000 DEFAULT 20 (Sets minimum intron length. Default: 20)
# PARAMETER OPTIONAL max.intron.length: "Maximum intron length" TYPE INTEGER FROM 1 TO 1000000 DEFAULT 500000 (Sets maximum intron length. Default: 500000)
# PARAMETER OPTIONAL library.type: "Library type" TYPE [fr-unstranded: fr-unstranded, fr-firststrand: fr-firststrand, fr-secondstrand: fr-secondstrand] DEFAULT fr-unstranded (Which library type to use. For directional\/strand specific library prepartion methods, choose fr-firststrand or fr-secondstrand depending on the preparation method: if the first read \(read1\) maps to the opposite, non-coding strand, choose fr-firststrand. If the first read maps to the coding strand, choose fr-secondstrand. For example for Illumina TruSeq Stranded sample prep, choose fr-firstsrand.)

# AO 30.5.2017 First version

# INPUT OPTIONAL reads1.txt: "List of read 1 files" TYPE GENERIC
# INPUT OPTIONAL reads2.txt: "List of read 2 files" TYPE GENERIC
# INPUT OPTIONAL genes.gtf: "Optional GTF file" TYPE GENERIC
# OUTPUT OPTIONAL hisat.bam.bai
# OUTPUT OPTIONAL junctions.bed
# OUTPUT OPTIONAL hisat-summary.txt
# PARAMETER OPTIONAL no.discordant: "Report only concordant alignments" TYPE [yes, no] DEFAULT yes (Report only concordant mappings.) 
# TODO: comment library.type away, unstranded is the default on hisat2
# PARAMETER OPTIONAL mate.inner.distance: "Expected inner distance between mate pairs" TYPE INTEGER DEFAULT 200 (Expected mean inner distance between mate pairs. For example, if your fragment size is 300 bp and read length is 50 bp, the inner distance is 200.)
# TODO: commonet inner.distance away? maxins is by default 500
# PARAMETER OPTIONAL use.gtf: "Use internal annotation GTF" TYPE [yes, no] DEFAULT yes (If this option is selected, TopHat will extract the transcript sequences and use Bowtie to align reads to this virtual transcriptome first. Only the reads that do not fully map to the transcriptome will then be mapped on the genome. The reads that did map on the transcriptome will be converted to genomic mappings (spliced as needed\) and merged with the novel mappings and junctions in the final TopHat output. If user provides a GTF file it is used instead of the internal annotation.)
# PARAMETER OPTIONAL no.novel.juncs: "When GTF file is used, ignore novel junctions" TYPE [yes, no] DEFAULT no (Only look for reads across junctions indicated in the supplied GTF file.)
# PARAMETER OPTIONAL mate.std.dev: "Standard deviation for the inner distances between mate pairs" TYPE INTEGER DEFAULT 20 (The standard deviation for the distribution on inner distances between mate pairs. The default is 20bp.)
# PARAMETER OPTIONAL max.multihits: "How many hits is a read allowed to have" TYPE INTEGER FROM 1 TO 1000000 DEFAULT 20 (Instructs TopHat to allow up to this many alignments to the reference for a given read.)
# PARAMETER OPTIONAL mismatches: "Number of mismatches allowed in final alignment" TYPE INTEGER FROM 0 TO 5 DEFAULT 2 (Final read alignments having more than this many mismatches are discarded.)
# PARAMETER OPTIONAL min.anchor.length: "Minimum anchor length" TYPE INTEGER FROM 3 TO 1000 DEFAULT 8 (TopHat2 will report junctions spanned by reads with at least this many bases on each side of the junction. Note that individual spliced alignments may span a junction with fewer than this many bases on one side. However, every junction involved in spliced alignments is supported by at least one read with this many bases on each side.)
# PARAMETER OPTIONAL splice.mismatches: "Maximum number of mismatches allowed in the anchor" TYPE INTEGER FROM 0 TO 2 DEFAULT 0 (The maximum number of mismatches that may appear in the anchor region of a spliced alignment.)
# PARAMETER OPTIONAL no.mixed: "Report only paired alignments" TYPE [yes, no] DEFAULT yes (Only report read alignments if both reads in a pair can be mapped.)

## Source required functions
source(file.path(chipster.common.path, "zip-utils.R"))
source(file.path(chipster.common.path, "tool-utils.R"))

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

## Parse parameters and store them into hisat.parameters
hisat.parameters <- ""
# Reads to align
hisat.parameters <-paste(hisat.parameters, "-U", "read.fq" )
# Quality score
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
hisat.parameters <- paste(hisat.parameters, "-x", organism)
# Set environment variable that specifies what index to use, HISAT2 requires this
Sys.setenv(HISAT2_INDEXES = "/opt/chipster/tools/genomes/indexes/hisat2")

## Set hard coded parameters
# Threads that hisat uses
hisat.parameters <- paste(hisat.parameters, "-p" , chipster.threads.max)
# Name of the output file
hisat.parameters <- paste(hisat.parameters,  "-S", "hisat.sam")
# Forward errors to hisat.log
hisat.parameters <- paste(hisat.parameters,  "2>> hisat.log")

#Print the HISAT2_INDEXES into debug
debugPrint("")
debugPrint("HISAT2_INDEXES:")
debugPrint("$HISAT2_INDEXES")




# Threads



# Print parameters into log
debugPrint(toString(hisat.parameters))

# Get input name
input.names <- read.table("chipster-inputs.tsv", header=F, sep="\t")

# check out if the file is compressed and if so unzip it
unzipInputs(input.names)

# setting up HISAT binaries (and paths)
hisat.binary <- file.path(chipster.tools.path, "hisat2", "hisat2")
samtools.binary <- file.path(chipster.tools.path, "samtools", "samtools")



## Run HISAT
#Parameters:-p Launch NTHREADS parallel search threads (default: 1). Threads will run on separate processors/cores and synchronize when parsing reads 
#				and outputting alignments. Searching for alignments is highly parallel, and speedup is close to linear. Increasing -p increases HISAT2's
#				memory footprint. E.g. when aligning to a human genome index, increasing -p from 1 to 8 increases the memory footprint by a few hundred
#				megabytes. This option is only available if bowtie is linked with the pthreads library (i.e. if BOWTIE_PTHREADS=0 is not specified at 
#				build time).
#           -q Reads (specified with <m1>, <m2>, <s>) are FASTQ files. FASTQ files usually have extension .fq or .fastq. FASTQ is the default format.
#				See also: --solexa-quals and --int-quals.
#           -x The basename of the index for the reference genome. The basename is the name of any of the index files up to but not including the 
#				final .1.ht2 / etc. hisat2 looks for the specified index first in the current directory, then in the directory specified in the 
#				HISAT2_INDEXES environment variable.
# 			-S File to write SAM alignments to. By default, alignments are written to the "standard out" or "stdout" filehandle (i.e. the console).
# 			-U Comma-separated list of files containing unpaired reads to be aligned, 
#				e.g. lane1.fq,lane2.fq,lane3.fq,lane4.fq. Reads may be a mix of different lengths.
#				If - is specified, hisat2 gets the reads from the "standard in" or "stdin" filehandle.
#			2>> Redirects error messages from stderr to spesified location.
# Note a single ' at the beginning, it allows us to use special characters like >
command <- paste("bash -c '", hisat.binary)

# Add the parameters
command <- paste(command, hisat.parameters)

# Close the command with a ', because there is a opening ' also
command <- paste(command,"'") 
# Print the command to the hisat.log file
system(paste("echo ",command ," >> hisat.log"))

# Run command
system(command)

## Run samtools
# Convert SAM file into BAM file and index bam file
# Parameters:
# 	-b Output in the BAM format.
#	-S Ignored for compatibility with previous samtools versions. Previously this option was required if input was in SAM format, but now the correct
#		format is automatically detected by examining the first few characters of input.
# Convert SAM to BAM
system(paste(samtools.binary, "view -bS hisat.sam > hisat.bam"))
# index bam
system(paste(samtools.binary, "sort hisat.bam"))
# Unset environmet variable
Sys.unsetenv("HISAT2_INDEX")

# Append the debug.log into hisat.log
system("cat debug.log >> hisat.log")

#EOF

