# TOOL hisat2-paired-end-with-index-building.R: "HISAT2 for paired end reads and own genome" (Aligns paired end RNA-seq reads to a user-supplied genome. You need to supply the genome in either as a FASTA format sequence or as a tar package with a HISAT2 index.)
# INPUT reads{...}.fq: "Reads to align" TYPE GENERIC
# INPUT OPTIONAL reads1.txt: "List of read 1 files" TYPE GENERIC
# INPUT OPTIONAL reads2.txt: "List of read 2 files" TYPE GENERIC
# INPUT OPTIONAL genome.txt: "Genome to align against \(fasta or HISAT2 index tar\)" TYPE GENERIC
# OUTPUT OPTIONAL hisat.bam
# OUTPUT OPTIONAL hisat.bam.bai
# OUTPUT OPTIONAL hisat.log
# OUTPUT OPTIONAL hisat2_index.tar
# PARAMETER rna.strandness: "RNA-strandness" TYPE [unstranded, FR, RF] DEFAULT unstranded (Specify strand-specific information. FR means read 1 is always on the same strand as the gene. RF means read 2 is always on the same strand as the gene. The default is unstranded.)
# PARAMETER quality.format: "Base quality encoding used" TYPE [sanger: "Sanger - Phred+33", phred64: "Phred+64"] DEFAULT sanger (Quality encoding used in the fastq file.)
# PARAMETER OPTIONAL max.multihits: "How many hits to report per read" TYPE INTEGER FROM 1 TO 1000000 DEFAULT 5 (Instructs HISAT2 to report up to this many alignments to the reference for a given read.)
# PARAMETER OPTIONAL min.intron.length: "Minimum intron length" TYPE INTEGER FROM 1 TO 1000 DEFAULT 20 (Sets minimum intron length. Default: 20)
# PARAMETER OPTIONAL max.intron.length: "Maximum intron length" TYPE INTEGER FROM 1 TO 1000000 DEFAULT 500000 (Sets maximum intron length. Default: 500000)
# PARAMETER OPTIONAL no.softclip: "Disallow soft clipping" TYPE [nosoft: "No soft-clipping", yessoft: "Use soft-clipping"] DEFAULT yessoft (Is soft-clipping used. By default HISAT2 may soft-clip reads near their 5' and 3' ends.)
# PARAMETER OPTIONAL dta: "Require long anchor lengths for subsequent assembly" TYPE [nodta: "Don't require", yesdta: "Require"] DEFAULT nodta (With this option, HISAT2 requires longer anchor lengths for de novo discovery of splice sites. This leads to fewer alignments with short-anchors, which helps transcript assemblers improve significantly in computation and memory usage.)
# RUNTIME R-4.1.1

# AO 30.5.2017 First version
# EK 18.10.2017 Polishing

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

# setting up HISAT binaries (and paths)
hisat.binary <- file.path(chipster.tools.path,"hisat2","hisat2")
samtools.binary <- c(file.path(chipster.tools.path, "samtools", "bin", "samtools"))
hisat.index.binary <- file.path(chipster.tools.path,"hisat2","hisat2-build")

# Get input name
input.names <- read.table("chipster-inputs.tsv",header = FALSE,sep = "\t")
input.display.names <- read_input_definitions()

# check out if the file is compressed and if so unzip it
unzipInputs(input.names)

if (fileOk("genome.txt")) {
  genome.filetype <- system("file -b genome.txt | cut -d ' ' -f2",intern = TRUE)
  hg_ifn <- ("")
  echo.command <- paste("echo Host genome file type",genome.filetype," > hisat.log")
  runExternal(echo.command)
  new_index_created <- ("no")
  # case 1. Ready calculated indexes in tar format
  if (genome.filetype == "tar") {
    runExternal("echo Extracting tar formatted gemome index file >> hisat.log")
    # Untar. Folders are flattened
    runExternal("tar xf genome.txt --xform='s#^.+/##x' 2>> hisat.log")
    # Check index base name
    if (file.exists(Sys.glob("*.1.ht2"))) {
      f <- list.files(getwd(),pattern = "\\.1.ht2$")
      hisat2.genome <- substr(f[1],1,nchar(f[1]) - 6)
    } else if (file.exists(Sys.glob("*.1.ht2l"))) {
      f <- list.files(getwd(),pattern = "\\.1.ht2l$")
      hisat2.genome <- substr(f[1],1,nchar(f[1]) - 7)
    } else {
      stop("CHIPSTER-NOTE: The .tar package does not seem to contain a valid Hisat2 index.")
    }
    # case 2. Fasta file
  } else {
    # Do indexing
    runExternal("echo Indexing the genome... >> hisat.log")
    runExternal("echo >> hisat.log")
    hisat2.genome <- strip_name(input.display.names$genome.txt)
    index.command <- paste("bash -c '",hisat.index.binary,"-p",chipster.threads.max,"genome.txt",hisat2.genome,"2>>hisat.log","'")
    runExternal(paste("echo ",index.command," >> hisat.log"))
    runExternal(index.command)
    echo.command <- paste("echo Internal genome name:",hisat2.genome," >> hisat.log")
    runExternal(echo.command)
    # Make tar package from the index
    if (file.exists(Sys.glob("*.1.ht*"))) {
      system("tar cf hisat2_index.tar *.ht*")
    }
  }
} else {
  stop("CHIPSTER-NOTE: Genome file not provided. Please make sure genome file is selcted and assigned to the correct input.")
}
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
} else if (file.exists("reads002.fq") && !file.exists("reads003.fq")) {
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
} else if (rna.strandness == "fr-secondstrand") {
  hisat.parameters <- paste(hisat.parameters,"--rna-strandness RF")
}
# Organism
# -x declares the basename of the index for reference genome, here we are using own genome
hisat.parameters <- paste(hisat.parameters,"-x",hisat2.genome)
# Set environment variable that defines where indexes locate, HISAT2 requires this
# Sys.setenv(HISAT2_INDEXES = "/opt/chipster/tools/genomes/indexes/hisat2")
# We have created own genome -> no need fot environment variable
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

## Set parameters that are not mutable via Chipster
# Threads that hisat uses
hisat.parameters <- paste(hisat.parameters,"-p",chipster.threads.max)
# Name of the output file
hisat.parameters <- paste(hisat.parameters,"-S","hisat.sam")
# Forward errors to hisat.log
hisat.parameters <- paste(hisat.parameters,"2>> hisat.log")
# Suppress SAM records for reads that failed to align
hisat.parameters <- paste(hisat.parameters,"--no-unal")

#Print the HISAT2_INDEXES into debug
debugPrint("")
debugPrint("HISAT2_INDEXES:")
debugPrint("$HISAT2_INDEXES")

# Print parameters into log
debugPrint("HISAT PARAMETERS")
debugPrint(toString(hisat.parameters))

## Run HISAT
# Note a single ' at the beginning, it allows us to use special characters like >
command <- paste("bash -c '",hisat.binary)

# Add the parameters
command <- paste(command,hisat.parameters)

# Close the command with a ', because there is a opening ' also
command <- paste(command,"'")
# Print the command to the hisat.log file
debugPrint(command)

# Run command
runExternal(command)

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
runExternal(samtools.view.command)
# Index bam, this produces a "hisat.sorted.bam" file
samtools.sort.command <- paste(samtools.binary,"sort hisat.tmp.bam -o hisat.sorted.bam")
debugPrint(samtools.sort.command)
runExternal(samtools.sort.command)

# Do not return empty BAM files
if (fileOk("hisat.sorted.bam",minsize = 100)) {
  # Rename result files
  runExternal("mv hisat.sorted.bam hisat.bam")
  # Change file names in BAM header to display names
  displayNamesToBAM("hisat.bam")
  # Index BAM
  runExternal(paste(samtools.binary,"index hisat.bam > hisat.bam.bai"))
}


# Unset environmet variable
Sys.unsetenv("HISAT2_INDEX")

# Append the debug.log into hisat.log
runExternal("cat debug.log >> hisat.log")

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
outputnames <- matrix(NA,nrow = 3,ncol = 2)
outputnames[1,] <- c("hisat.bam",paste(basename,".bam",sep = ""))
outputnames[2,] <- c("hisat.bam.bai",paste(basename,".bam.bai",sep = ""))
outputnames[3,] <- c("hisat2_index.tar",paste(hisat2.genome,".hisat2.tar",sep = ""))

# Write output definitions file
write_output_definitions(outputnames)

#EOF
