# TOOL bowtie-paired-end-with-index-building.R: "Bowtie for paired end reads and own genome" (This tool uses Bowtie to align paired-end reads to a provided reference genome. You need to supply the paired-end reads as two FASTQ files containing the reads in the same order. You also need to provide the reference genome as input in FASTA format.)
# INPUT reads1.fq: "No 1 mate reads" TYPE FASTQ
# INPUT reads2.fq: "No 2 mate reads" TYPE FASTQ
# INPUT genome.txt: "Genome to align against" TYPE GENERIC
# OUTPUT bowtie.bam
# OUTPUT bowtie.bam.bai
# OUTPUT bowtie.log
# OUTPUT OPTIONAL unaligned_1.fq
# OUTPUT OPTIONAL unaligned_2.fq
# OUTPUT OPTIONAL multireads_1.fq
# OUTPUT OPTIONAL multireads_2.fq
# PARAMETER max.mismatches: "Number of mismatches allowed" TYPE [0, 1, 2, 3] DEFAULT 2 (How many mismatches are the alignments allowed to have?)
# PARAMETER limit.to.seed: "Consider mismatches only in the seed region" TYPE [yes, no] DEFAULT no (Should the mismatch limit be applied only to the left, good quality part of the read? You can define the length of this seed region with the next parameter.)
# PARAMETER seed: "Length of the seed region" TYPE INTEGER FROM 5 TO 50 DEFAULT 28 (If you have chosen to apply the mismatch limit only to the left, good quality part of the read, how many bases should be considered? The minimum length of seed region is 5.)
# PARAMETER multiread: "How many places is a read allowed to align to" TYPE [1, 2, 1000000: "no limit"] DEFAULT 1000000 (If you want to have alignments only for uniquely mapping reads, select 1.)
# PARAMETER OPTIONAL min.insert.size: "Minimum insert size" TYPE INTEGER FROM 0 TO 1000 DEFAULT 0 (The minimum insert size for valid paired-end alignments. E.g. if 60 is specified and a paired-end alignment consists of two 20-bp alignments in the appropriate orientation with a 20-bp gap between them, that alignment is considered valid.)
# PARAMETER OPTIONAL max.insert.size: "Maximum insert size" TYPE INTEGER FROM 50 TO 1500 DEFAULT 250 (The maximum insert size for valid paired-end alignments. E.g. if 100 is specified and a paired-end alignment consists of two 20-bp alignments in the proper orientation with a 60-bp gap between them, that alignment is considered valid.)
# PARAMETER OPTIONAL orientation: "Upstream-downstream mate orientation" TYPE [fr: "mate1 upstream of reverse complement of mate2 or vice versa", rf: "upstream mate1 reverse-complemented and mate2 forward-oriented"] DEFAULT fr (The upstream-downstream mate orientations for a valid paired-end alignment against the forward reference strand.)
# PARAMETER OPTIONAL quality: "Allowed total of mismatch qualities" TYPE INTEGER FROM 10 TO 100 DEFAULT 70 (What is the maximum permitted total of quality values of ALL mismatch positions throughout the read (not just in the seed region\)? Note that this parameter is taken into account only if you have chosen to apply the mismatch limit to the seed region.)
# PARAMETER OPTIONAL quality.format: "Quality value format used" TYPE [solexa1_3: "Illumina GA v1.3-1.5", sanger: Sanger] DEFAULT sanger (Note that this parameter is taken into account only if you chose to apply the mismatch limit to the seed region. Are the quality values in the Sanger format (ASCII characters equal to the Phred quality plus 33\) or in the Illumina Genome Analyzer Pipeline v1.3 or later format (ASCII characters equal to the Phred quality plus 64\)? Please see the manual for details.)
# PARAMETER OPTIONAL alignment.no: "How many valid alignments are reported per read" TYPE [1, 2, 3] DEFAULT 1 (If there are several, equally good alignments, how many should be reported?)
# PARAMETER OPTIONAL multiread.file: "Put multireads to a separate file" TYPE [yes, no] DEFAULT no (If you chose not to have alignments for reads which map to multiple positions, would you like to store these reads to a separate fastq file?)
# PARAMETER OPTIONAL unaligned.file: "Put unaligned reads to a separate file" TYPE [yes, no] DEFAULT no (Would you like to store unaligned reads to a new fastq file? Note that also multireads will be added to this file, unless you asked them to be put to a separate file.)
# RUNTIME R-4.1.1


# EK 12.7.2011
# AMS 19.6.2012 Added unzipping
# EK 1.11.2012 fixed genome parameter and SAM output
# AMS 11.11.2013 Added thread support

source(file.path(chipster.common.lib.path, "tool-utils.R"))

# check out if the file is compressed and if so unzip it
source(file.path(chipster.common.lib.path, "zip-utils.R"))
unzipIfGZipFile("reads1.fq")
unzipIfGZipFile("reads2.fq")
unzipIfGZipFile("genome.txt")


# bowtie
bowtie.binary <- c(file.path(chipster.tools.path, "bowtie", "bowtie"))
version <- system(paste(bowtie.binary, "--version | head -1 | cut -d ' ' -f 3"), intern = TRUE)
documentVersion("Bowtie", version)
bowtie.index.binary <- c(file.path(chipster.tools.path, "bowtie", "bowtie-build"))

genome.filetype <- system("file -b genome.txt | cut -d ' ' -f2", intern = TRUE)
print(paste("Host genome file type", genome.filetype))

# case 1. Ready calculated indexes in tar format
if (genome.filetype == "tar") {
  print("Extracting tar formatted gemome index file")
  runExternal("tar xf genome.txt")
  # Check index base name
  if (length(Sys.glob("*.1.ebwt")) != 0) {
    f <- list.files(getwd(), pattern = "\\.1.ebwt$")
    bowtie.genome <- substr(f[1], 1, nchar(f[1]) - 7)
  } else {
    stop("CHIPSTER-NOTE: The .tar package does not seem to contain a valid Bowtie index.")
  }

  runExternal("ls -l")
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

  # bowtie index building
  print("Indexing the genome...")
  bowtie.genome <- "genome"

  command.indexing <- paste(bowtie.index.binary, "-f", "genome.txt", bowtie.genome)
  runExternal(command.indexing)
  runExternal("ls -l")
}
print(paste("Internal genome name:", bowtie.genome))
command.start <- paste("bash -c '", bowtie.binary)

# common parameters
common.parameters <- paste("-p", chipster.threads.max, "-S", "-q", "-m", multiread, "-k", alignment.no, "-I", min.insert.size, "-X", max.insert.size)

# mode specific parameters
quality.parameter <- ifelse(quality.format == "solexa1_3", "--solexa1.3-quals", "")
orientation.parameter <- ifelse(orientation == "rf", "--rf", "")
n.mode.parameters <- paste("-n", max.mismatches, "-l", seed, "-e", quality, quality.parameter)
v.mode.parameters <- paste("-v", max.mismatches)
mode.parameters <- ifelse(limit.to.seed == "yes", n.mode.parameters, v.mode.parameters)

# output parameters
unaligned.output <- ifelse(unaligned.file == "yes", "--un unaligned.fq", "")
multiread.output <- ifelse(multiread.file == "yes", "--max multireads.fq", "")
output.parameters <- paste(unaligned.output, multiread.output)

# command ending
command.end <- paste("-x", bowtie.genome, "-1 reads1.fq -2 reads2.fq 1> alignment.sam 2> bowtie.log'")

# run bowtie
bowtie.command <- paste(command.start, common.parameters, quality.parameter, orientation.parameter, mode.parameters, output.parameters, command.end)
# stop(paste('CHIPSTER-NOTE: ', bowtie.command))
runExternal(bowtie.command)


# samtools binary
samtools.binary <- c(file.path(chipster.tools.path, "samtools", "bin", "samtools"))
version <- system(paste(samtools.binary, "--version | head -1 | cut -d ' ' -f 2"), intern = TRUE)
documentVersion("SAMtools", version)


# convert sam to bam
runExternal(paste(samtools.binary, "view -bS -q 1 alignment.sam -o alignment.bam"))

# sort bam
runExternal(paste(samtools.binary, "sort alignment.bam -o alignment.sorted.bam"))

# index bam
runExternal(paste(samtools.binary, "index alignment.sorted.bam"))

# rename result files
runExternal("mv alignment.sorted.bam bowtie.bam")
runExternal("mv alignment.sorted.bam.bai bowtie.bam.bai")

# Handle output names
#
source(file.path(chipster.common.lib.path, "tool-utils.R"))

# read input names
inputnames <- read_input_definitions()

# Determine base name
base1 <- strip_name(inputnames$reads1.fq)
base2 <- strip_name(inputnames$reads2.fq)
basename <- paired_name(base1, base2)

# Make a matrix of output names
outputnames <- matrix(NA, nrow = 6, ncol = 2)
outputnames[1, ] <- c("bowtie.bam", paste(basename, ".bam", sep = ""))
outputnames[2, ] <- c("bowtie.bam.bai", paste(basename, ".bam.bai", sep = ""))
outputnames[3, ] <- c("unaligned_1.fq", paste(base1, "_unaligned.fq", sep = ""))
outputnames[4, ] <- c("unaligned_2.fq", paste(base2, "_unaligned.fq", sep = ""))
outputnames[5, ] <- c("multireads_1.fq", paste(base1, "_multireads.fq", sep = ""))
outputnames[6, ] <- c("multireads_2.fq", paste(base2, "_multireads.fq", sep = ""))

# Write output definitions file
write_output_definitions(outputnames)
