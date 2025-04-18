# TOOL bwa-mem-paired-end-with-index-building.R: "BWA MEM for paired-end reads and own genome" (Aligns reads to genomes using the BWA MEM algorithm. Results are sorted and indexed BAM files.
# Note that this BWA MEM tool requires that you have imported the reference genome to Chipster in fasta format. If you would like to align paired-end reads against publicly available genomes, please use the tool \"BWA MEM for paired-end reads\".)
# INPUT reads1.txt: "Reads to align" TYPE FASTQ
# INPUT reads2.txt: "Reads to align" TYPE FASTQ
# INPUT genome.txt: "Reference genome" TYPE GENERIC
# OUTPUT bwa.bam
# OUTPUT bwa.log
# OUTPUT OPTIONAL bwa.bam.bai
# PARAMETER OPTIONAL index.file: "Create index file" TYPE [index_file: "Create index file", no_index: "No index file"] DEFAULT no_index (Creates index file for BAM. By default no index file.)
# PARAMETER mode: "Data source" TYPE [ normal: " Illumina, 454, IonTorrent reads longer than 70 base pairs", pacbio: "PacBio subreads"] DEFAULT normal (Defining the type of reads will instruct the tool to use a predefined set of parameters optimized for that read type.)
# RUNTIME R-4.1.1

# KM 11.11.2014

# check out if the file is compressed and if so unzip it
source(file.path(chipster.common.lib.path, "zip-utils.R"))
unzipIfGZipFile("reads1.txt")
unzipIfGZipFile("reads2.txt")
unzipIfGZipFile("genome.txt")

# bwa
bwa.binary <- file.path(chipster.tools.path, "bwa", "bwa mem")
bwa.index.binary <- file.path(chipster.module.path, "shell", "check_bwa_index.sh")
command.start <- paste("bash -c '", bwa.binary)

# Do indexing
print("Indexing the genome...")
runExternal("echo Indexing the genome... > bwa.log")
check.command <- paste(bwa.index.binary, "genome.txt| tail -1 ")

# genome.dir <- system(check.command, intern = TRUE)
# bwa.genome <- file.path( genome.dir , "genome.txt")
bwa.genome <- system(check.command, intern = TRUE)

mode.parameters <- ifelse(mode == "pacbio", "-x pacbio", "")

# command ending
command.end <- paste(bwa.genome, "reads1.txt reads2.txt 1> alignment.sam 2>> bwa.log'")

# run bwa alignment
bwa.command <- paste(command.start, mode.parameters, command.end)

echo.command <- paste("echo '", bwa.binary, mode.parameters, bwa.genome, "reads.txt ' > bwa.log")
# stop(paste('CHIPSTER-NOTE: ', bwa.command))
runExternal(echo.command)
runExternal(bwa.command)

# samtools binary
samtools.binary <- c(file.path(chipster.tools.path, "samtools", "bin", "samtools"))

# convert sam to bam
runExternal(paste(samtools.binary, "view -bS alignment.sam -o alignment.bam"))

# sort bam
runExternal(paste(samtools.binary, "sort alignment.bam -o alignment.sorted.bam"))

# index bam
runExternal(paste(samtools.binary, "index alignment.sorted.bam"))

# rename result files
runExternal("mv alignment.sorted.bam bwa.bam")
if (index.file == "index_file") {
  runExternal("mv alignment.sorted.bam.bai bwa.bam.bai")
}

# Handle output names
#
source(file.path(chipster.common.lib.path, "tool-utils.R"))

# read input names
inputnames <- read_input_definitions()

# Determine base name
base1 <- strip_name(inputnames$reads1.txt)
base2 <- strip_name(inputnames$reads2.txt)
basename <- paired_name(base1, base2)

# Make a matrix of output names
outputnames <- matrix(NA, nrow = 2, ncol = 2)
outputnames[1, ] <- c("bwa.bam", paste(basename, ".bam", sep = ""))
outputnames[2, ] <- c("bwa.bam.bai", paste(basename, ".bam.bai", sep = ""))

# Write output definitions file
write_output_definitions(outputnames)
