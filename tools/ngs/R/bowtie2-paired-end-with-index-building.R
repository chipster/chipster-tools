# TOOL bowtie2-paired-end-with-index-building.R: "Bowtie2 for paired end reads and own genome" (This tool uses Bowtie2 to align paired-end reads to a provided reference genome. Paired-end reads can either be in FASTA or FASTQ format, but both reads files need to be in the same format. You also need to provide the reference genome as input either in FASTA format or as a tar package with a Bowtie2 index.)
# INPUT reads{...}.fq: "Reads" TYPE FASTQ
# INPUT OPTIONAL reads1.txt: "List of read 1 files" TYPE TEXT_STRICT
# INPUT OPTIONAL reads2.txt: "List of read 2 files" TYPE TEXT_STRICT
# INPUT genome.txt: "Genome to align against" TYPE GENERIC
# OUTPUT bowtie2.bam
# OUTPUT bowtie2.log
# OUTPUT OPTIONAL bowtie2.bam.bai
# OUTPUT OPTIONAL failed_1.fq
# OUTPUT OPTIONAL failed_2.fq
# OUTPUT OPTIONAL bowtie2_index.tar
# PARAMETER OPTIONAL index.file: "Create index file" TYPE [index_file: "Create index file", no_index: "No index file"] DEFAULT no_index (Creates index file for BAM. By default no index file.)
# PARAMETER strategy: "Alignment strategy to use" TYPE [--very-fast: "Very fast", --fast: "Fast", --sensitive: "Sensitive", --very-sensitive: "Very sensitive", --very-fast-local: "Very fast local", -fast-local: "Fast local", --sensitive-local: "Sensitive local", --very-sensitive-local: "Very sensitive local"] DEFAULT --sensitive (The alignment strategy to be used. Bowtie2 can map the reads using end-to-end or local alignments. When local alignment is used, Bowtie2 might "trim" or "clip" some read characters from one or both ends of the alignment if doing so maximizes the alignment score. Bowtie2 uses heuristics for mapping the reads to the reference genome. Several Bowtie2 parameters affect simultaneously both to the sensitivity and to computing time. In Chipster you can choose the sensitivity level from a set of pre-defined parameter combinations that allow you to tune the balance between the computing time and mapping sensitivity.)
# PARAMETER quality.format: "Quality value format used" TYPE [--phred33: "Sanger - Phred+33", --phred64: "Illumina GA v1.3-1.5 - Phred+64", --ignore-quals: "Fixed 30 for all"] DEFAULT --phred33 (Quality scale used in the fastq-file.)
# PARAMETER alignment.no: "How many valid alignments are reported per read" TYPE [0: "Best based on the mapping quality", 1: "1", 2: "2", 3: "3", 4: "4", 5: "5", 6: "All alignments"] DEFAULT 0 (By default, Bowtie2 reports only the best alignment of the read (based on the mapping quality\). Optionally, if there are several, equally good alignments, you can choose how many of them should be reported?)
# PARAMETER OPTIONAL discordant.file: "Put reads that did not align concordantly to a separate file" TYPE [yes, no] DEFAULT no (Write paired-end reads that fail to align concordantly to a fastq file. This includes reads that aligned discordantly, reads whose mate failed to align and unaligned reads.)
# PARAMETER OPTIONAL ma: "Match bonus" TYPE INTEGER FROM 0 TO 10 DEFAULT 2 (Match bonus for a match in local alignment. Default value 2)
# PARAMETER OPTIONAL mp: "Maximum penalty for mismatch" TYPE INTEGER FROM 0 TO 20 DEFAULT 6 (Maximum penalty for mismatch; lower quality = lower penalty. Default value 6)
# PARAMETER OPTIONAL np: "Penalty for non-ACGTs"  TYPE INTEGER FROM 0 TO 20 DEFAULT 1 ( Sets penalty for positions where the read, reference, or both, contain an ambiguous character such as N. Default: 1.)
# PARAMETER OPTIONAL rdg.open: "Gap opening penalty for the reads" TYPE INTEGER FROM 0 TO 20 DEFAULT 5 (Gap opening penalty for the reads. Default value: 5. )
# PARAMETER OPTIONAL rdg.ext: "Gap extension penalty for the reads" TYPE INTEGER FROM 0 TO 20 DEFAULT 3 (Gap extension penalty for the reads. Default value: 3. )
# PARAMETER OPTIONAL rfg.open: "Gap opening penalty for the reference" TYPE INTEGER FROM 0 TO 20 DEFAULT 5 (Gap opening penalty for the reference. Default value: 5. )
# PARAMETER OPTIONAL rfg.ext: "Gap extension penalty for the reference" TYPE INTEGER FROM 0 TO 20 DEFAULT 3 (Gap extension penalty for the reference. Default value: 3. )
# PARAMETER OPTIONAL minins: "Minimum insert length" TYPE INTEGER FROM 0 TO 2000 DEFAULT 0 (Minimum insert length between the mate pairs. Default value: 0)
# PARAMETER OPTIONAL maxins: "Maximum insert length" TYPE INTEGER FROM 0 TO 4000 DEFAULT 500 (Maximum insert length between the mate pairs. Default value: 500)
# PARAMETER OPTIONAL pair.order: "Order of mates to align" TYPE [--fr: "Forward/reverse", --rf: "Reverse/Forward", --ff: "Forward/forward"] DEFAULT --fr (The orientation of the mate pairs. Default: forward/revrse)
# PARAMETER OPTIONAL no.mixed: "Suppress unpaired alignments" TYPE [yes, no] DEFAULT no (By default, when bowtie2 cannot find a concordant or discordant alignment for a pair, it then tries to find alignments for the individual mates. This option disables that behavior.)
# PARAMETER OPTIONAL no.discordant: "Suppress discordant alignments" TYPE [yes, no] DEFAULT no (By default, bowtie2 looks for discordant alignments if it cannot find any concordant alignments. A discordant alignment is an alignment where both mates align uniquely, but that does not satisfy the paired-end constraints. This option disables that behavior)
# PARAMETER OPTIONAL no.dovetail: "Not concordant when mates extend past each other" TYPE [yes, no] DEFAULT no (If the mates "dovetail", that is if one mate alignment extends past the beginning of the other such that the wrong mate begins upstream, consider that to be concordant. Default: mates cannot dovetail in a concordant alignment. )
# PARAMETER OPTIONAL no.contain: "Not concordant when one mate alignment contains other" TYPE [yes, no] DEFAULT no (If one mate alignment contains the other, consider that to be non-concordant. Default: a mate can contain the other in a concordant alignment.)
# PARAMETER OPTIONAL no.overlap: "Not concordant when mates overlap at all"  TYPE [yes, no] DEFAULT no (If one mate alignment overlaps the other at all, consider that to be non-concordant. Default: mates can overlap in a concordant alignment.)
# RUNTIME R-4.1.1

# KM 10-01.2012
# EK 8.5.2013 replaced samtools -q 1 with Bowtie --no-unal to remove unaligned reads from BAM
# AMS 11.11.2013 Added thread support
# When updating Bowtie2 to 2.2.x, remember to change mp parameter

# PARAMETER OPTIONAL unaligned.file: "Put unaligned reads to a separate file" TYPE [yes, no] DEFAULT no (Would you like to store unaligned reads to a new fastq file? Note that also multireads will be added to this file, unless you asked them to be put to a separate file.)

source(file.path(chipster.common.lib.path, "tool-utils.R"))
source(file.path(chipster.common.lib.path, "bam-utils.R"))
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
  runExternal("echo Extracting tar formatted gemome index file >> bowtie2.log")
  runExternal("tar -tf genome.txt >> bowtie2.log")
  check.command <- paste(bowtie2.index.binary, "genome.txt | tail -1 ")
  bowtie2.genome <- system(check.command, intern = TRUE)
  runExternal("ls -l >> bowtie2.log")
  # case 2. Fasta file
} else {
  # Do indexing
  # check sequence file type
  emboss.path <- file.path(chipster.tools.path, "emboss-20.04", "bin")
  options(scipen = 999)
  inputfile.to.check <- ("genome.txt")
  sfcheck.binary <- file.path(chipster.module.path, "../misc/shell/sfcheck.sh")
  sfcheck.command <- paste(sfcheck.binary, emboss.path, inputfile.to.check)
  str.filetype <- system(sfcheck.command, intern = TRUE)

  if (str.filetype == "Not an EMBOSS compatible sequence file") {
    stop("CHIPSTER-NOTE: Your reference genome is not a sequence file that is compatible with the tool you try to use")
  }


  print("Indexing the genome...")
  runExternal("echo Indexing the genome... >> bowtie2.log")
  check.command <- paste(bowtie2.index.binary, "genome.txt -tar | tail -1 ")
  bowtie2.genome <- system(check.command, intern = TRUE)
  runExternal("ls -l >> bowtie2.log")
  new_index_created <- ("yes")
}
echo.command <- paste("echo Internal genome name:", bowtie2.genome, " >> bowtie2.log")
runExternal(echo.command)

command.start <- paste("bash -c '", bowtie.binary)
rdg.value <- paste(rdg.open, rdg.ext, sep = ",")
rfg.value <- paste(rfg.open, rfg.ext, sep = ",")

parameters <- paste(strategy, "--mp", mp, "--np", np, "--rdg", rdg.value, "--rfg", rfg.value, "--minins", minins, "--maxins", maxins, pair.order, quality.format, "--no-unal", "-p", chipster.threads.max)

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

if (no.mixed == "yes") {
  parameters <- paste(parameters, "--no-mixed")
}

if (no.discordant == "yes") {
  parameters <- paste(parameters, "--no-discordant")
}

if (no.dovetail == "yes") {
  parameters <- paste(parameters, "--no-dovetail")
}

if (no.contain == "yes") {
  parameters <- paste(parameters, "--no-contain")
}

if (no.overlap == "yes") {
  parameters <- paste(parameters, "--no-overlap")
}

if (discordant.file == "yes") {
  parameters <- paste(parameters, "--un-conc failed")
}

# Check if reads are in FASTA format
# emboss.path <- file.path(chipster.tools.path,"emboss","bin")
# sfcheck.binary <- file.path(chipster.module.path,"../misc/shell/sfcheck.sh")
# sfcheck.command <- paste(sfcheck.binary,emboss.path,"reads001.fq")
# str.filetype <- system(sfcheck.command,intern = TRUE)
# if (str.filetype == "fasta") {
#   parameters <- paste(parameters,"-f")
# }

# Input files
if (file.exists("reads1.txt") && file.exists("reads2.txt")) {
  # Case: list files exist
  reads1.list <- make_input_list("reads1.txt")
  reads2.list <- make_input_list("reads2.txt")
  if (identical(intersect(reads1.list, reads2.list), character(0))) {
    reads1 <- paste(reads1.list, sep = "", collapse = ",")
    reads2 <- paste(reads2.list, sep = "", collapse = ",")
  } else {
    stop(paste("CHIPSTER-NOTE: ", "One or more files is listed in both lists."))
  }
} else if (file.exists("reads002.fq") && !file.exists("reads003.fq")) {
  # Case: no list file, but only two fastq inputs
  in.sorted <- input.names[order(input.names[, 2]), ]
  reads <- grep("reads", in.sorted[, 1], value = TRUE)
  reads1 <- reads[1]
  reads2 <- reads[2]
} else {
  # Case: no list files, more than two fastq inputs
  stop(paste("CHIPSTER-NOTE: ", "List file is missing. You need to provide a list of read files for both directions."))
}

# output parameters
# output.parameters <- paste(unaligned.output, multiread.output)
# stop(paste('CHIPSTER-NOTE: ', parameters))
# command ending
command.end <- paste("-x", bowtie2.genome, "-1", reads1, "-2", reads2, "1> alignment.sam 2>> bowtie2.log'")

# run bowtie
bowtie.command <- paste(command.start, parameters, command.end)
# stop(paste('CHIPSTER-NOTE: ', bowtie.command))

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

# Change file named in BAM header to display names
displayNamesToBAM("alignment.bam")

# sort bam
runExternal(paste(samtools.binary, "sort alignment.bam -o alignment.sorted.bam"))

# index bam
runExternal(paste(samtools.binary, "index alignment.sorted.bam"))

# rename result files according to the index parameter
runExternal("mv alignment.sorted.bam bowtie2.bam")
if (index.file == "index_file") {
  runExternal("mv alignment.sorted.bam.bai bowtie2.bam.bai")
}

# if (unaligned.file== "yes"){
#  system("mv unaligned.1 unaligned_1.fq")
#  system("mv unaligned.2 unaligned_2.fq")
# }

if (discordant.file == "yes") {
  runExternal("mv failed.1 failed_1.fq")
  runExternal("mv failed.2 failed_2.fq")
}

# Substitute display names to log for clarity
displayNamesToFile("bowtie2.log")

# Handle output names
#
# read input names
inputnames <- read_input_definitions()

# Determine base name
name1 <- unlist(strsplit(reads1, ","))
base1 <- strip_name(inputnames[[name1[1]]])

name2 <- unlist(strsplit(reads2, ","))
base2 <- strip_name(inputnames[[name2[1]]])

basename <- paired_name(base1, base2)

# Make a matrix of output names
outputnames <- matrix(NA, nrow = 5, ncol = 2)
outputnames[1, ] <- c("bowtie2.bam", paste(basename, ".bam", sep = ""))
outputnames[2, ] <- c("bowtie2.bam.bai", paste(basename, ".bam.bai", sep = ""))
outputnames[3, ] <- c("failed_1.fq", paste(base1, "_failed.fq", sep = ""))
outputnames[4, ] <- c("failed_2.fq", paste(base2, "_failed.fq", sep = ""))
if (new_index_created == "yes") {
  hg_ifn <- strip_name(inputnames$genome.txt)
  outputnames[5, ] <- c("bowtie2_index.tar", paste(hg_ifn, "_bowtie2_index.tar", sep = ""))
}

# Write output definitions file
write_output_definitions(outputnames)
