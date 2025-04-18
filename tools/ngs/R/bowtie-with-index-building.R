# TOOL bowtie-with-index-building.R: "Bowtie for single end reads and own genome" (This tool uses Bowtie to align single-end reads to a provided reference genome. You need to supply the single-end reads in FASTQ format. You also need to provide the reference genome as input in FASTA format.)
# INPUT reads.txt: "Reads to align" TYPE FASTQ
# INPUT genome.txt: "Genome to align against" TYPE GENERIC
# OUTPUT bowtie.bam
# OUTPUT bowtie.bam.bai
# OUTPUT bowtie.log
# OUTPUT OPTIONAL unaligned-reads.fastq
# OUTPUT OPTIONAL multireads.fastq
# PARAMETER max.mismatches: "Number of mismatches allowed" TYPE [0, 1, 2, 3] DEFAULT 2 (How many mismatches are the alignments allowed to have?)
# PARAMETER limit.to.seed: "Consider mismatches only in the seed region" TYPE [yes, no] DEFAULT no (Should the mismatch limit be applied only to the left, good quality part of the read? You can define the length of this seed region with the next parameter.)
# PARAMETER seed: "Length of the seed region" TYPE INTEGER FROM 5 TO 50 DEFAULT 28 (If you have chosen to apply the mismatch limit only to the left, good quality part of the read, how many bases should be considered? The minimum length of seed region is 5.)
# PARAMETER quality: "Allowed total of mismatch qualities" TYPE INTEGER FROM 10 TO 100 DEFAULT 70 (What is the maximum permitted total of quality values of ALL mismatch positions throughout the read (not just in the seed region\)? Note that this parameter is taken into account only if you have chosen to apply the mismatch limit to the seed region.)
# PARAMETER quality.format: "Quality value format used" TYPE [solexa1_3: "Illumina GA v1.3-1.5", sanger: Sanger] DEFAULT sanger (Note that this parameter is taken into account only if you chose to apply the mismatch limit to the seed region. Are the quality values in the Sanger format (ASCII characters equal to the Phred quality plus 33\) or in the Illumina Genome Analyzer Pipeline v1.3 or later format (ASCII characters equal to the Phred quality plus 64\)? Please see the manual for details.)
# PARAMETER OPTIONAL multiread: "How many best category hits is the read allowed to have" TYPE [1, 2, 1000000: "no limit"] DEFAULT 1000000 (If you select 1, only those reads which have a unique best category hit will be reported. Note that the read can still have other, lower category hits even if they are not reported. The hit categories are based on the number of mismatches.)
# PARAMETER OPTIONAL alignment.no: "How many valid alignments are reported per read" TYPE [1, 2, 3] DEFAULT 1 (If there are several, equally good alignments, how many should be reported?)
# PARAMETER OPTIONAL multiread.file: "Put multireads to a separate file" TYPE [yes, no] DEFAULT no (If you chose not to have alignments for reads which map to multiple positions, would you like to store these reads to a separate fastq file?)
# PARAMETER OPTIONAL unaligned.file: "Put unaligned reads to a separate file" TYPE [yes, no] DEFAULT no (Would you like to store unaligned reads to a new fastq file? Note that also multireads will be added to this file, unless you asked them to be put to a separate file.)
# RUNTIME R-4.1.1

# EK 20.6.2011
# AMS 19.6.2012 Added unzipping
# AMS 11.11.2013 Added thread support

source(file.path(chipster.common.lib.path, "tool-utils.R"))

# check out if the file is compressed and if so unzip it
source(file.path(chipster.common.lib.path, "zip-utils.R"))
unzipIfGZipFile("reads.txt")
unzipIfGZipFile("genome.txt")




# bowtie
bowtie.binary <- c(file.path(chipster.tools.path, "bowtie", "bowtie"))
version <- system(paste(bowtie.binary, "--version | head -1 | cut -d ' ' -f 3"), intern = TRUE)
documentVersion("Bowtie", version)
bowtie.index.binary <- c(file.path(chipster.tools.path, "bowtie", "bowtie-build"))

genome.filetype <- system("file -b genome.txt | cut -d ' ' -f2", intern = TRUE)
print(paste("echo Host genome file type", genome.filetype))

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
common.parameters <- paste("-p", chipster.threads.max, "-q --best -S --strata", "-m", multiread, "-k", alignment.no)

# mode specific parameters
quality.parameter <- ifelse(quality.format == "solexa1_3", "--solexa1.3-quals", "")
n.mode.parameters <- paste("-n", max.mismatches, "-l", seed, "-e", quality, quality.parameter)
v.mode.parameters <- paste("-v", max.mismatches)
mode.parameters <- ifelse(limit.to.seed == "yes", n.mode.parameters, v.mode.parameters)

# output parameters
unaligned.output <- ifelse(unaligned.file == "yes", "--un unaligned-reads.fastq", "")
multiread.output <- ifelse(multiread.file == "yes", "--max multireads.fastq", "")
output.parameters <- paste(unaligned.output, multiread.output)

# command ending
command.end <- paste("-x", bowtie.genome, "reads.txt 1> alignment.sam 2> bowtie.log'")

# run bowtie
bowtie.command <- paste(command.start, common.parameters, quality.parameter, mode.parameters, output.parameters, command.end)
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
basename <- strip_name(inputnames$reads.txt)

# Make a matrix of output names
outputnames <- matrix(NA, nrow = 4, ncol = 2)
outputnames[1, ] <- c("bowtie.bam", paste(basename, ".bam", sep = ""))
outputnames[2, ] <- c("bowtie.bam.bai", paste(basename, ".bam.bai", sep = ""))
outputnames[3, ] <- c("unaligned-reads.fastq", paste(basename, "_unaligned.fq", sep = ""))
outputnames[4, ] <- c("multireads.fastq", paste(basename, "_multireads.fq", sep = ""))

# Write output definitions file
write_output_definitions(outputnames)
