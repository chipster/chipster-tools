# TOOL hisat2-paired-end-test.R: "HISAT2 for paired end reads - small footprint test" (Aligns paired end RNA-seq reads to a reference genome.)
# INPUT reads{...}.fq.gz: "Reads to align" TYPE GENERIC
# INPUT OPTIONAL reads1.txt: "List of read 1 files" TYPE GENERIC
# INPUT OPTIONAL reads2.txt: "List of read 2 files" TYPE GENERIC
# OUTPUT OPTIONAL hisat.bam
# OUTPUT OPTIONAL hisat.bam.bai
# OUTPUT OPTIONAL hisat.log
# PARAMETER organism: "Genome" TYPE ["FILES genomes/indexes/hisat2 .fa"] DEFAULT "SYMLINK_TARGET genomes/indexes/hisat2/default .fa" (Genome or transcriptome that you would like to align your reads against.)
# PARAMETER rna.strandness: "RNA-strandness" TYPE [unstranded, FR, RF] DEFAULT unstranded (Specify strand-specific information. FR means read 1 is always on the same strand as the gene. RF means read 2 is always on the same strand as the gene. The default is unstranded.)
# PARAMETER quality.format: "Base quality encoding used" TYPE [sanger: "Sanger - Phred+33", phred64: "Phred+64"] DEFAULT sanger (Quality encoding used in the fastq file.)
# PARAMETER OPTIONAL max.multihits: "How many hits to report per read" TYPE INTEGER FROM 1 TO 1000000 DEFAULT 5 (Instructs HISAT2 to report up to this many alignments to the reference for a given read.)
# PARAMETER OPTIONAL min.intron.length: "Minimum intron length" TYPE INTEGER FROM 1 TO 1000 DEFAULT 20 (Sets minimum intron length. Default: 20)
# PARAMETER OPTIONAL max.intron.length: "Maximum intron length" TYPE INTEGER FROM 1 TO 1000000 DEFAULT 500000 (Sets maximum intron length. Default: 500000)
# PARAMETER OPTIONAL no.softclip: "Disallow soft-clipping" TYPE [nosoft: "No soft-clipping", yessoft: "Use soft-clipping"] DEFAULT yessoft (Is soft-cliping used. By default HISAT2 may soft-clip reads near their 5' and 3' ends.)
# PARAMETER OPTIONAL dta: "Require long anchor lengths for subsequent assembly" TYPE [nodta: "Don't require", yesdta: "Require"] DEFAULT nodta (With this option, HISAT2 requires longer anchor lengths for de novo discovery of splice sites. This leads to fewer alignments with short-anchors, which helps transcript assemblers improve significantly in computation and memory usage.)
# PARAMETER OPTIONAL dta.cufflinks: "Tailor alignments for Cufflinks" TYPE [yes, no] DEFAULT no (Report alignments tailored specifically for Cufflinks transcript assembly. With this option, HISAT2 looks for novel splice sites with three signals, GT-AG, GC-AG, AT-AC, but all user-provided splice sites are used irrespective of their signals. HISAT2 produces an optional XS field for every spliced alignment.)

# AO 30.05.2017 First version
# EK 18.10.2017 Polishing
# EK 29.11.2021 Add --dta-cufflinks option

## Source required functions
source(file.path(chipster.common.path,"zip-utils.R"))
source(file.path(chipster.common.path,"tool-utils.R"))
source(file.path(chipster.common.path,"bam-utils.R"))

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
    system(paste("echo ",command,">> debug.log"))
  }
}

## Options
# Prefer fixed representation over exponential
options(scipen = 10)
# Debug mode, change debug to TRUE or FALSE, depending do you want debug prints or not
debug <- FALSE
debugPrint("")
debugPrint("DEBUG MODE IS ON")

# Get input name
input.names <- read.table("chipster-inputs.tsv",header = FALSE,sep = "\t")

# check out if the file is compressed and if so unzip it
#unzipInputs(input.names)

## Parse parameters and store them into hisat.parameters
hisat.parameters <- ""
# Reads to align
# Parse the read names from input files
if (file.exists("reads1.txt") && file.exists("reads2.txt")) {
  # Case: list files exist
  reads1.list <- make_input_list("reads1.txt")
  reads2.list <- make_input_list("reads2.txt")
  if (identical(intersect(reads1.list,reads2.list),character(0))) {
    reads1 <- paste(reads1.list,sep = "",collapse = ",")
    reads2 <- paste(reads2.list,sep = "",collapse = ",")
  } else {
    stop(paste('CHIPSTER-NOTE: ',"One or more files is listed in both lists."))
  }
} else if (file.exists("reads002.fq.gz") && !file.exists("reads003.fq.gz")) {
  # Case: no list file, but only two fastq inputs
  in.sorted <- input.names[order(input.names[,2]),]
  reads <- grep("reads",in.sorted[,1],value = TRUE)
  reads1 <- reads[1]
  reads2 <- reads[2]
} else {
  # Case: no list files, more than two fastq inputs
  stop(paste('CHIPSTER-NOTE: ',"List file is missing. You need to provide a list of read files for both directions."))
}
hisat.parameters <- paste(hisat.parameters,"-1",reads1,"-2",reads2)
# Quality score format
if (quality.format == "phred64") {
  hisat.parameters <- paste(hisat.parameters,"--phred64")
} else {
  hisat.parameters <- paste(hisat.parameters,"--phred33")
}
# Intron length, defaluts are 20 and 500 000
hisat.parameters <- paste(hisat.parameters,"--min-intronlen",min.intron.length)
hisat.parameters <- paste(hisat.parameters,"--max-intronlen",max.intron.length)
# Specify strand-specific information: the default is unstranded
if (rna.strandness == "FR") {
  hisat.parameters <- paste(hisat.parameters,"--rna-strandness FR")
} else if (rna.strandness == "RF") {
  hisat.parameters <- paste(hisat.parameters,"--rna-strandness RF")
}
# Organism
# -x declares the basename of the index for reference genome
hisat.parameters <- paste(hisat.parameters,"-x",organism)
# Set environment variable that defines where indexes locate, HISAT2 requires this
Sys.setenv(HISAT2_INDEXES = "/opt/chipster/tools/genomes/indexes/hisat2")
# Known splice sites
if (file.exists("splicesites.txt")) {
  hisat.parameters <- paste(hisat.parameters,"--known-splicesite-infile","splicesites.txt")
}
# How many hits is a read allowed to have
hisat.parameters <- paste(hisat.parameters,"-k",max.multihits)
# Allow soft-clipping, by default soft-clipping is used
if (no.softclip == "nosoft") {
  hisat.parameters <- paste(hisat.parameters,"--no-softclip")
}
# Is longer anchor lengths required
if (dta == "yesdta") {
  hisat.parameters <- paste(hisat.parameters,"--dta")
}

# Should alignments be trailored for Cufflinks
if (dta.cufflinks == "yes") {
	hisat.parameters <- paste(hisat.parameters,"--dta-cufflinks")
}

## Set parameters that are not mutable via Chipster
# Threads that hisat uses
hisat.parameters <- paste(hisat.parameters,"-p",chipster.threads.max)
# Name of the output file
#hisat.parameters <- paste(hisat.parameters,"-S","hisat.sam")
# Forward errors to hisat.log
#hisat.parameters <- paste(hisat.parameters,"2>> hisat.log")
# Suppress SAM records for reads that failed to align
hisat.parameters <- paste(hisat.parameters,"--no-unal")

#Print the HISAT2_INDEXES into debug
debugPrint("")
debugPrint("HISAT2_INDEXES:")
debugPrint("$HISAT2_INDEXES")

# Print parameters into log
debugPrint("HISAT PARAMETERS")
debugPrint(toString(hisat.parameters))

# setting up HISAT binaries (and paths)
hisat.binary <- file.path(chipster.tools.path,"hisat2","hisat2")
samtools.binary <- file.path(chipster.tools.path,"samtools","samtools")
hisat.index <- file.path(chipster.tools.path, "genomes", "indexes", "hisat2", organism)

## Run HISAT
# Note a single ' at the beginning, it allows us to use special characters like >
command <- paste(hisat.binary, "-x", hisat.index)

# Add the parameters
command <- paste(command,hisat.parameters, "2> hisat.log |", samtools.binary, "view -bS - > hisat.tmp.bam")

# Close the command with a ', because there is a opening ' also
#command <- paste(command,"'")
documentCommand(command)
# Print the command to the hisat.log file
debugPrint(command)

# Run command
system(command)

## Run samtools
# Convert SAM file into BAM file and index bam file
# Parameters:
#   -b Output in the BAM format.
#   -S Ignored for compatibility with previous samtools versions. Previously this option was required if input was in SAM format, but now the correct
#    format is automatically detected by examining the first few characters of input.
# Convert SAM to BAM
debugPrint("")
debugPrint("SAMTOOLS")
samtools.view.command <- paste(samtools.binary,"view -bS hisat.sam > hisat.tmp.bam")
debugPrint(samtools.view.command)
#system(samtools.view.command)
# Index bam, this produces a "hisat.sorted.bam" file
samtools.sort.command <- paste(samtools.binary,"sort hisat.tmp.bam hisat.sorted")
debugPrint(samtools.sort.command)
system(samtools.sort.command)

# Do not return empty BAM files
if (fileOk("hisat.sorted.bam",minsize = 100)) {
  # Rename result files
  system("mv hisat.sorted.bam hisat.bam")
  # Change file names in BAM header to display names
  displayNamesToBAM("hisat.bam")
  # Index BAM
  system(paste(samtools.binary,"index hisat.bam > hisat.bam.bai"))
}


# Unset environmet variable
Sys.unsetenv("HISAT2_INDEX")

# Append the debug.log into hisat.log
#system("cat debug.log >> hisat.log")

# Substitute display names to log for clarity
displayNamesToFile("hisat.log")

# Handle output names
#

# read input names
inputnames <- read_input_definitions()

# Determine base name
name1 <- unlist(strsplit(reads1,","))
base1 <- strip_name(inputnames[[name1[1]]])

name2 <- unlist(strsplit(reads2,","))
base2 <- strip_name(inputnames[[name2[1]]])

basename <- paired_name(base1,base2)

# Make a matrix of output names
outputnames <- matrix(NA,nrow = 2,ncol = 2)
outputnames[1,] <- c("hisat.bam",paste(basename,".bam",sep = ""))
outputnames[2,] <- c("hisat.bam.bai",paste(basename,".bam.bai",sep = ""))

# Write output definitions file
write_output_definitions(outputnames)

system("ls | sort -r >> hisat.log")
#EOF

