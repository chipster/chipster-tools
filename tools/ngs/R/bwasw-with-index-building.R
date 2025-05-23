# TOOL bwasw-with-index-building.R: "BWA-SW for single end reads and own genome" (BWA aligns reads to genomes and transcriptomes. This BWA version uses the BWA-SW algorithm that is designed for long, over 300pb, low-quality reads. Results are sorted and indexed bam files.
# Note that this BWA tool requires that you have imported the reference genome to Chipster in fasta format. If you would like to align reads against publicly available genomes, please use the tool \"BWA-SW for single end reads\".)
# INPUT reads.txt: "Reads to align" TYPE FASTQ
# INPUT genome.txt: "Reference genome" TYPE GENERIC
# OUTPUT bwa.bam
# OUTPUT bwa.bam.bai
# OUTPUT bwa.log
# PARAMETER OPTIONAL match.score: "Score of a match" TYPE INTEGER DEFAULT 1 (Score of a matching base. Corresponds to the command line parameter -a)
# PARAMETER OPTIONAL gap.opening: "Gap opening penalty " TYPE INTEGER DEFAULT 11 (Gap opening penalty. Corresponds to the command line parameter -q)
# PARAMETER OPTIONAL gap.extension: "Gap extension penalty " TYPE INTEGER DEFAULT 4 (Gap extension penalty. Corresponds to the command line parameter -r)
# PARAMETER OPTIONAL mismatch.penalty: "Mismatch penalty" TYPE INTEGER DEFAULT 3 (Mismatch penalty. Corresponds to the command line parameter -b)
# PARAMETER OPTIONAL band.width: "Band width in the banded alignment"  TYPE INTEGER DEFAULT 33 (Band width in the banded alignment. Corresponds to the command line parameter -w .)
# PARAMETER OPTIONAL min.score: "Minimum score threshold divided by the given value"  TYPE INTEGER DEFAULT 37 (Minimum score threshold divided by the given value. Corresponds to the command line parameter -T.)
# PARAMETER OPTIONAL threshold.coeff: "Coefficient for threshold adjustment" TYPE DECIMAL DEFAULT 5.5 (Coefficient for threshold adjustment according to query length. Given an l-long query, the threshold for a hit to be retained is a*max\{T,c*log\(l\)\}. Corresponds to the command line parameter -c. )
# PARAMETER OPTIONAL z.best: "Z-best heuristics" TYPE INTEGER DEFAULT 1 (Z-best heuristics. Higher -z increases accuracy at the cost of speed. Corresponds to the command line parameter -z.)
# PARAMETER OPTIONAL sa.interval: "Maximum SA interval size" TYPE INTEGER DEFAULT 3 (Maximum SA interval size for initiating a seed. Higher value increases accuracy at the cost of speed. Corresponds to the command line parameter -s.)
# PARAMETER OPTIONAL min.support: "Reverse alignment limit" TYPE INTEGER DEFAULT 5 (Minimum number of seeds supporting the resultant alignment to skip reverse alignment. Corresponds to the command line parameter -N)
# PARAMETER OPTIONAL alignment.no: "How many valid alignments are reported per read" TYPE  INTEGER DEFAULT 3 (Maximum number of alignments to report. Corresponds to the command line parameter bwa samse -n )
# RUNTIME R-4.1.1

# KM 24.8.2011
# AMS 19.6.2012 Added unzipping
# AMS 11.11.2013 Added thread support

# check out if the file is compressed and if so unzip it
source(file.path(chipster.common.lib.path, "zip-utils.R"))
unzipIfGZipFile("reads.txt")
unzipIfGZipFile("genome.txt")

# bwa
bwa.binary <- file.path(chipster.tools.path, "bwa", "bwa")
bwa.index.binary <- file.path(chipster.module.path, "shell", "check_bwa_index.sh")
command.start <- paste("bash -c '", bwa.binary)

# Do indexing
print("Indexing the genome...")
runExternal("echo Indexing the genome... > bwa.log")
check.command <- paste(bwa.index.binary, "genome.txt| tail -1 ")
# genome.dir <- system(check.command, intern = TRUE)
# bwa.genome <- file.path( genome.dir , "genome.txt")
bwa.genome <- system(check.command, intern = TRUE)

# algorithm parameters
mode.parameters <- paste("bwasw", "-t", chipster.threads.max, "-b", mismatch.penalty, "-q", gap.opening, "-r", gap.extension, "-a", match.score, "-w", band.width, "-T", min.score, "-c", threshold.coeff, "-z", z.best, "-s", sa.interval, "-N", min.support)

# command ending
command.end <- paste(bwa.genome, "reads.txt 1> alignment.sai 2>> bwa.log'")

# run bwa alignment
runExternal("echo Running the alignment with command: >> bwa.log")
bwa.command <- paste(command.start, mode.parameters, command.end)
echo.command <- paste("echo '", bwa.command, "'  >> bwa.log")
runExternal(echo.command)
# stop(paste('CHIPSTER-NOTE: ', bwa.command))
runExternal(bwa.command)

# sai to sam conversion
runExternal("echo Running sai to sam conversion with command: >> bwa.log")
samse.parameters <- paste("samse -n", alignment.no)
samse.end <- paste(bwa.genome, "-f alignment.sam alignment.sai reads.txt >> bwa.log'")
samse.command <- paste(command.start, samse.parameters, samse.end)
echo.command <- paste("echo '", samse.command)
runExternal(echo.command)
echo.command <- paste("echo '", samse.end, " >> bwa.log")
runExternal(echo.command)
runExternal(samse.command)


# samtools binary
samtools.binary <- c(file.path(chipster.tools.path, "samtools", "bin", "samtools"))

# convert sam to bam
runExternal("echo Converting sam to bam format. >> bwa.log")
runExternal(paste(samtools.binary, "view -bS alignment.sam -o alignment.bam"))

# sort bam
runExternal(paste(samtools.binary, "sort alignment.bam -o alignment.sorted.bam"))

# index bam
runExternal(paste(samtools.binary, "index alignment.sorted.bam"))

# rename result files
runExternal("mv alignment.sorted.bam bwa.bam")
runExternal("mv alignment.sorted.bam.bai bwa.bam.bai")
