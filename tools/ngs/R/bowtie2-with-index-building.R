# TOOL bowtie2-with-index-building.R: "Bowtie2 for single end reads and own genome" (This tool uses Bowtie2 to align single-end reads to a provided reference genome. Single-end reads can either be in FASTA or FASTQ format. You also need to provide the reference genome as input either in FASTA format or as a tar package with a Bowtie2 index.)
# INPUT reads{...}.fq: "Reads to align" TYPE FASTQ
# INPUT genome.txt: "Genome to align against" TYPE GENERIC
# OUTPUT OPTIONAL bowtie2.bam
# OUTPUT OPTIONAL bowtie2.bam.bai
# OUTPUT bowtie2.log
# OUTPUT OPTIONAL unaligned_1.fq
# OUTPUT OPTIONAL bowtie2_index.tar
# PARAMETER OPTIONAL index.file: "Create index file" TYPE [index_file: "Create index file", no_index: "No index file"] DEFAULT no_index (Creates index file for BAM. By default no index file.)
# PARAMETER strategy: "Alignment strategy to use" TYPE [--very-fast: "Very fast", --fast: "Fast", --sensitive: "Sensitive", --very-sensitive: "Very sensitive", --very-fast-local: "Very fast local", -fast-local: "Fast local", --sensitive-local: "Sensitive local", --very-sensitive-local: "Very sensitive local"] DEFAULT --sensitive (The alignment strategy to be used. Bowtie2 can map the reads using end-to-end or local alignments. When local alignment is used, Bowtie2 might "trim" or "clip" some read characters from one or both ends of the alignment if doing so maximizes the alignment score. Bowtie2 uses heuristics for mapping the reads to the reference genome. Several Bowtie2 parameters affect simultaneously both to the sensitivity and to computing time. In Chipster you can choose the sensitivity level with a set of pre-defined parameter combinations that allow you to tune the balance between the computing time and mapping sensitivity.)
# PARAMETER quality.format: "Quality value format used" TYPE [--phred33: "Sanger - Phred+33", --phred64: "Illumina GA v1.3-1.5 - Phred+64", --ignore-quals: "Fixed 30 for all"] DEFAULT --phred33 (Quality scale used in the fastq-file.)
# PARAMETER alignment.no: "How many valid alignments are reported per read" TYPE [0: "Best based on the mapping quality", 1: "1", 2: "2", 3: "3", 4: "4", 5: "5", 6: "All alignments"] DEFAULT 0 (By default, Bowtie2 reports only the best alignment of the read (based on the mapping quality\). Optionally, if there are several, equally good alignments, you can choose how many of them should be reported?)
# PARAMETER OPTIONAL unaligned.file: "Put unaligned reads to a separate file" TYPE [yes, no] DEFAULT no (Store unaligned reads to a new fastq file.)
# PARAMETER OPTIONAL ma: "Match bonus" TYPE INTEGER FROM 0 TO 10 DEFAULT 2 (Match bonus for a match in local alignment. Default value 2)
# PARAMETER OPTIONAL mp: "Maximum penalty for mismatch" TYPE INTEGER FROM 0 TO 20 DEFAULT 6 (Maximum penalty for mismatch; lower quality = lower penalty. Default value 6)
# PARAMETER OPTIONAL np: "Penalty for non-ACGTs"  TYPE INTEGER FROM 0 TO 20 DEFAULT 1 ( Sets penalty for positions where the read, reference, or both, contain an ambiguous character such as N. Default: 1.)
# PARAMETER OPTIONAL rdg.open: "Gap opening penalty for the reads" TYPE INTEGER FROM 0 TO 20 DEFAULT 5 (Gap opening penalty for the reads. Default value: 5. )
# PARAMETER OPTIONAL rdg.ext: "Gap extension penalty for the reads" TYPE INTEGER FROM 0 TO 20 DEFAULT 3 (Gap extension penalty for the reads. Default value: 3. )
# PARAMETER OPTIONAL rfg.open: "Gap opening penalty for the reference" TYPE INTEGER FROM 0 TO 20 DEFAULT 5 (Gap opening penalty for the reference. Default value: 5. )
# PARAMETER OPTIONAL rfg.ext: "Gap extension penalty for the reference" TYPE INTEGER FROM 0 TO 20 DEFAULT 3 (Gap extension penalty for the reference. Default value: 3. )
# RUNTIME R-4.1.1

# KM 23.10.2012
# EK 8.5.2013 replaced samtools -q 1 with Bowtie --no-unal to remove unaligned reads from BAM
# AMS 11.11.2013 Added thread support
# When updating Bowtie2 to 2.2.x, remember to change mp parameter

source(file.path(chipster.common.lib.path, "bam-utils.R"))
source(file.path(chipster.common.lib.path, "tool-utils.R"))
source(file.path(chipster.common.lib.path, "zip-utils.R"))

genomeFile <- "genome.txt"
if (!file.exists(genomeFile)) {
  stop("CHIPSTER-NOTE: Reference genome is not defined! Please assign value for parameter: Genome to align against")
}

# check out if the file is compressed and if so unzip it
input.names <- read.table("chipster-inputs.tsv", header = FALSE, sep = "\t")
for (i in 1:nrow(input.names)) {
  unzipIfGZipFile(input.names[i, 1])
}

# bowtie2
bowtie.binary <- c(file.path(chipster.tools.path, "bowtie2", "bowtie2"))
version <- system(paste(bowtie.binary, "--version | head -1 | cut -d ' ' -f 3"), intern = TRUE)
documentVersion("Bowtie", version)
bowtie2.index.binary <- file.path(chipster.module.path, "shell", "check_bowtie2_index.sh")


genome.filetype <- system("file -b genome.txt | cut -d ' ' -f2", intern = TRUE)
hg_ifn <- ("")
echo.command <- paste("echo Host genome file type", genome.filetype, " > bowtie2.log")
runExternal(echo.command)

new_index_created <- ("no")
# case 1. Ready calculated indexes in tar format
if (genome.filetype == "tar") {
  runExternal("echo Extarting tar formatted gemome index file >> bowtie2.log")
  runExternal("tar -tf genome.txt >> bowtie2.log")
  check.command <- paste(bowtie2.index.binary, "genome.txt | tail -1 ")
  bowtie2.genome <- system(check.command, intern = TRUE)
  # system("ls -l >> bowtie2.log")
  # case 2. Fasta file
} else {
  # check sequece file type
  emboss.path <- file.path(chipster.tools.path, "emboss-20.04", "bin")
  options(scipen = 999)
  inputfile.to.check <- ("genome.txt")
  sfcheck.binary <- file.path(chipster.module.path, "../misc/shell/sfcheck.sh")
  sfcheck.command <- paste(sfcheck.binary, emboss.path, inputfile.to.check)
  str.filetype <- system(sfcheck.command, intern = TRUE)
  if (str.filetype == "Not an EMBOSS compatible sequence file") {
    stop("CHIPSTER-NOTE: Your reference genome is not a sequence file that is compatible with the tool you try to use")
  }

  # Do indexing
  print("Indexing the genome...")
  runExternal("echo Indexing the genome... >> bowtie2.log")
  check.command <- paste(bowtie2.index.binary, "genome.txt -tar | tail -1 ")
  bowtie2.genome <- system(check.command, intern = TRUE)
  # genome.dir <- system(check.command, intern = TRUE)
  # bowtie2.genome <- file.path( genome.dir , "genome.txt")
  new_index_created <- ("yes")
}
# echo.command <- paste("echo Internal genome name:", bowtie2.genome, " >> bowtie2.log")
# system(echo.command)
# echo.command <- paste("echo current >> bowtie2.log; ls -l >> bowtie2.log; echo ",bowtie2.genome, " >> bowtie2.log; ls -l ", bowtie2.genome, "* >> bowtie2.log", sep="")
# system(echo.command)
# bowtie.genome <- c(file.path(chipster.tools.path, "bowtie2", "indexes" , genome))
command.start <- paste("bash -c '", bowtie.binary)
rdg.value <- paste(rdg.open, rdg.ext, sep = ",")
rfg.value <- paste(rfg.open, rfg.ext, sep = ",")

parameters <- paste(strategy, "--mp", mp, "--np", np, "--rdg", rdg.value, "--rfg", rfg.value, quality.format, "--no-unal", "-p", chipster.threads.max)

if (alignment.no > 0) {
  if (alignment.no == 6) {
    parameters <- paste(parameters, "-a")
  }
  if (alignment.no < 6) {
    parameters <- paste(parameters, "-k", alignment.no)
  }
}

# Local alignment specific parameters
if (strategy == "--very-fast-local" || strategy == "--fast-local" || strategy == "--sensitive-local" || strategy == "--very-sensitive-local") {
  parameters <- paste(parameters, "--local --ma", ma)
}

if (unaligned.file == "yes") {
  parameters <- paste(parameters, "--un unaligned")
}


# output parameters
# output.parameters <- paste(unaligned.output, multiread.output)
# stop(paste('CHIPSTER-NOTE: ', parameters))

# Check if reads are in FASTA format
# emboss.path <- file.path(chipster.tools.path,"emboss","bin")
# sfcheck.binary <- file.path(chipster.module.path,"../misc/shell/sfcheck.sh")
# sfcheck.command <- paste(sfcheck.binary,emboss.path,"reads001.fq")
# str.filetype <- system(sfcheck.command,intern = TRUE)
# if (str.filetype == "fasta") {
#  parameters <- paste(parameters,"-f")
# }

# Input fastq names
reads1 <- paste(grep("reads", input.names[, 1], value = TRUE), sep = "", collapse = ",")

# command ending
command.end <- paste("-x", bowtie2.genome, "-U", reads1, "1> alignment.sam 2>> bowtie2.log'")

# run bowtie
bowtie.command <- paste(command.start, parameters, command.end)
# stop(paste('CHIPSTER-NOTE: ', bowtie.command))
runExternal("echo Launching Bowtie2 alignment >> bowtie2.log")
echo.command <- paste("echo '", bowtie.command, "' >> bowtie2.log")
runExternal(echo.command)
runExternal(bowtie.command)

if (file.size("alignment.sam") < 1) {
  runExternal("cat bowtie2.log")
  stop("Bowtie2 failed! Check the tail of the ouput below for more information.")
}
# samtools binary
samtools.binary <- c(file.path(chipster.tools.path, "samtools", "bin", "samtools"))
version <- system(paste(samtools.binary, "--version | head -1 | cut -d ' ' -f 2"), intern = TRUE)
documentVersion("SAMtools", version)

# convert sam to bam
runExternal(paste(samtools.binary, "view -bS alignment.sam -o alignment.bam"))

# sort bam
runExternal(paste(samtools.binary, "sort alignment.bam -o alignment.sorted.bam"))

# index bam
runExternal(paste(samtools.binary, "index alignment.sorted.bam"))

# Substitute display names to BAM header for clarity
displayNamesToBAM("alignment.sorted.bam")

# rename result files according to the index parameter
runExternal("mv alignment.sorted.bam bowtie2.bam")
if (index.file == "index_file") {
  runExternal("mv alignment.sorted.bam.bai bowtie2.bam.bai")
}

if (unaligned.file == "yes") {
  runExternal("mv unaligned unaligned_1.fq")
}

# Substitute display names to log for clarity
displayNamesToFile("bowtie2.log")

# Handle output names
#

# read input names
inputnames <- read_input_definitions()

# Determine base name
basename <- strip_name(inputnames$reads001.fq)

# Make a matrix of output names
outputnames <- matrix(NA, nrow = 4, ncol = 2)
outputnames[1, ] <- c("bowtie2.bam", paste(basename, ".bam", sep = ""))
outputnames[2, ] <- c("bowtie2.bam.bai", paste(basename, ".bam.bai", sep = ""))
outputnames[3, ] <- c("unaligned_1.fq", paste(basename, "_unaligned.fq", sep = ""))
if (new_index_created == "yes") {
  hg_ifn <- strip_name(inputnames$genome.txt)
  outputnames[4, ] <- c("bowtie2_index.tar", paste(hg_ifn, "_bowtie2_index.tar", sep = ""))
}
# Write output definitions file
write_output_definitions(outputnames)
