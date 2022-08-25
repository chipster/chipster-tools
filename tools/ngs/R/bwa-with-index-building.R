# TOOL bwa-with-index-building.R: "BWA-backtrack for single end reads and own genome" (BWA-backtrack aligns paired end reads to genomes with BWA ALN algorithm. If more than one FASTQ file is provided, each file is first aligned separately, and the BAM files are then merged. Results are sorted and indexed bam files. 
# Note that this BWA tool requires that you have imported the reference genome to Chipster in fasta format. If you would like to align reads against publicly available genomes, please use the tool \"BWA for single end reads\".)
# INPUT reads{...}.fq: "Reads to align" TYPE GENERIC 
# INPUT genome.txt: "Reference genome" TYPE GENERIC
# OUTPUT bwa.bam 
# OUTPUT bwa.bam.bai 
# OUTPUT OPTIONAL bwa.log 
# OUTPUT OPTIONAL genome_bwa_index.tar
# PARAMETER seed.length: "Length of the seed region" TYPE INTEGER DEFAULT 32 (How many bases of the left, good quality part of the read should be used as the seed region. If the seed length is longer than the reads, the seeding will be disabled. Corresponds to the command line parameter -l.) 
# PARAMETER seed.edit: "Maximum number of of differences in the seed region" TYPE INTEGER DEFAULT 2 (Maximum number of differences such as mismatches or indels in the seed region. Corresponds to the command line parameter -k.)
# PARAMETER total.edit: "Maximum edit distance for the whole read" TYPE DECIMAL DEFAULT 0.04 ( Maximum edit distance if the value is more than one. If the value is between 1 and 0 then it defines the fraction of missing alignments given 2% uniform base error rate. In the latter case, the maximum edit distance is automatically chosen for different read lengths. Corresponds to the command line parameter -n.)
# PARAMETER quality.format: "Quality value format used" TYPE [solexa1_3: "Illumina GA v1.3-1.5", sanger: Sanger] DEFAULT sanger (Note that this parameter is taken into account only if you chose to apply the mismatch limit to the seed region. Are the quality values in the Sanger format (ASCII characters equal to the Phred quality plus 33\) or in the Illumina Genome Analyzer Pipeline v1.3 or later format (ASCII characters equal to the Phred quality plus 64\)? Please see the manual for details.)
# PARAMETER OPTIONAL num.gaps: "Maximum number of gaps" TYPE INTEGER DEFAULT 1 (Maximum number of gap openings for one read. Corresponds to the command line parameter -o.)
# PARAMETER OPTIONAL num.extensions: "Maximum number of gap extensions" TYPE INTEGER DEFAULT -1 (Maximum number of gap extensions, -1 for disabling long gaps. Corresponds to the command line parameter -e.)
# PARAMETER OPTIONAL gap.opening: "Gap opening penalty " TYPE INTEGER DEFAULT 11 (Gap opening penalty. Corresponds to the command line parameter -O.)
# PARAMETER OPTIONAL gap.extension: "Gap extension penalty " TYPE INTEGER DEFAULT 4 (Gap extension penalty. Corresponds to the command line parameter -E.)
# PARAMETER OPTIONAL mismatch.penalty: "Mismatch penalty threshold" TYPE INTEGER DEFAULT 4 (BWA will not search for suboptimal hits with a score lower than the alignment score minus this. Corresponds to the command line parameter -M.)
# PARAMETER OPTIONAL disallow.gaps: "Maximum occurrences for extending a long deletion"  TYPE INTEGER DEFAULT 10 (Maximum occurrences for extending a long deletion. Corresponds to the command line parameter -d.)
# PARAMETER OPTIONAL disallow.indel: "Disallow an indel within the given number of bp towards the ends"  TYPE INTEGER DEFAULT 5 (Do not put an indel within the defined value of bp towards the ends. Corresponds to the command line parameter -i.)
# PARAMETER OPTIONAL trim.threshold: "Quality trimming threshold" TYPE INTEGER DEFAULT 0 (Quality threshold for read trimming down to 35bp. Corresponds to the command line parameter -q.)
# PARAMETER OPTIONAL barcode.length: "Barcode length"  TYPE INTEGER DEFAULT 0 (Length of barcode starting from the 5 prime end. The barcode of each read will be trimmed before mapping. Corresponds to the command line parameter -B.)
# PARAMETER OPTIONAL alignment.no: "How many valid alignments are reported per read" TYPE  INTEGER DEFAULT 3 (Maximum number of alignments to report. Corresponds to the samse command line parameter -n.)
# RUNTIME R-4.1.1

# KM 24.8.2011
# AMS 19.6.2012 Added unzipping
# AMS 11.11.2013 Added thread support

source(file.path(chipster.common.path, "tool-utils.R"))
source(file.path(chipster.common.path, "bam-utils.R"))
source(file.path(chipster.common.path, "zip-utils.R"))

# check out if the file is compressed and if so unzip it
input.names <- read.table("chipster-inputs.tsv", header=F, sep="\t")
for (i in 1:nrow(input.names)) {
	unzipIfGZipFile(input.names[i,1])	
}

# read input names
inputnames <- read_input_definitions()

# bwa settings
bwa.binary <- file.path(chipster.tools.path, "bwa", "bwa")
bwa.index.binary <- file.path(chipster.module.path, "shell", "check_bwa_index.sh")
samtools.binary <- c(file.path(chipster.tools.path, "samtools", "bin", "samtools"))

command.start <- paste("bash -c '", bwa.binary)


# Do indexing
if (file.exists("genome.txt")){
	bwa.index.binary <- file.path(chipster.module.path, "shell", "check_bwa_index.sh")
	genome.filetype <- system("file -b genome.txt | cut -d ' ' -f2", intern = TRUE )
	
	# case 1. Ready calculated indexes in tar format
	if (genome.filetype == "tar"){
		check.command <- paste ( bwa.index.binary, "genome.txt| tail -1 ")
		bwa.genome <- system(check.command, intern = TRUE)		
		# case 2. Fasta file
	}else{
		check.command <- paste ( bwa.index.binary, "genome.txt -tar| tail -1 ")
		bwa.genome <- system(check.command, intern = TRUE)
		cp.command <- paste("cp ", bwa.genome, "_bwa_index.tar ./genome_bwa_index.tar ", sep ="")
		runExternal(cp.command)
		hg_ifn <- strip_name(inputnames$genome.txt)
		# Make a matrix of output names
		outputnames <- matrix(NA, nrow=1, ncol=2)
		outputnames[1,] <- c("genome_bwa_index.tar", paste(hg_ifn, "_bwa_index.tar", sep =""))
		# Write output definitions file
		write_output_definitions(outputnames)
	}
}else{
	stop(paste('CHIPSTER-NOTE: ', "Reference genome not found."))
} 




print("Indexing the genome...")
runExternal("echo Indexing the genome... > bwa.log")
check.command <- paste ( bwa.index.binary, "genome.txt| tail -1 ")
#genome.dir <- system(check.command, intern = TRUE)
#bwa.genome <- file.path( genome.dir , "genome.txt")
bwa.genome <- system(check.command, intern = TRUE)


# mode specific parameters
if (total.edit >= 1) {
	total.edit <- round(total.edit)
}
quality.parameter <- ifelse(quality.format == "solexa1_3", "-I", "")
mode.parameters <- paste("aln", "-t", chipster.threads.max, "-o", num.gaps, "-e", num.extensions, "-d", disallow.gaps, "-i" , disallow.indel , "-l" , seed.length , "-k" , seed.edit , "-O" , gap.opening , "-E" , gap.extension , "-q" , trim.threshold, "-B" , barcode.length, "-M" , mismatch.penalty , "-n" , total.edit , quality.parameter)


# command ending
command.end <- paste( bwa.genome , "reads.txt 1> alignment.sai 2>> bwa.log'")
runExternal("echo Running the alignment with command: >> bwa.log")

# Run BWA for each input
for (i in 1:nrow(input.names)) {
	# command ending
	sai.file <- paste(c(as.character(i), ".sai"), collapse="")
	sam.file <- paste(c(as.character(i), ".sam"), collapse="")
	bam.file <- paste(c(as.character(i), ".bam"), collapse="")
	command.end <- paste(bwa.genome, input.names[i,1], "1>", sai.file, "2>> bwa.log'")
	
	# run bwa alignment
	bwa.command <- paste(command.start, mode.parameters, command.end)

	documentCommand(bwa.command)
	
	runExternal(bwa.command)
	
	# sai to sam conversion
	samse.parameters <- paste("samse -n", alignment.no )
	samse.end <- paste(bwa.genome, sai.file, input.names[i,1], ">", sam.file, "2>>bwa.log'" )
	samse.command <- paste( command.start, samse.parameters , samse.end )
	runExternal(samse.command)
	
	# convert sam to bam
	runExternal(paste(samtools.binary, "view -bS", sam.file, "-o", bam.file))
}

# Join bam files
if (fileOk("2.bam")){
	# more than one bam exists, so join them
	runExternal("ls *.bam > bam.list")
	runExternal(paste(samtools.binary, "merge -b bam.list alignment.bam"))
}else{
	# only one bam, so just rename it
	runExternal("mv 1.bam alignment.bam")
}

# Change file named in BAM header to display names
displayNamesToBAM("alignment.bam")

# sort bam
runExternal(paste(samtools.binary, "sort alignment.bam -o alignment.sorted.bam"))

# index bam
runExternal(paste(samtools.binary, "index alignment.sorted.bam"))

# rename result files
runExternal("mv alignment.sorted.bam bwa.bam")
runExternal("mv alignment.sorted.bam.bai bwa.bam.bai")

# Substitute display names to log for clarity
displayNamesToFile("bwa.log")

## Handle output names

# read input names
inputnames <- read_input_definitions()

if (fileNotOk("reads002.fq")){
	# Determine base name
	basename <- strip_name(inputnames$reads001.fq)
}else{
	basename <- "bwa_multi"
}
# Make a matrix of output names
outputnames <- matrix(NA, nrow=2, ncol=2)
outputnames[1,] <- c("bwa.bam", paste(basename, ".bam", sep =""))
outputnames[2,] <- c("bwa.bam.bai", paste(basename, ".bam.bai", sep =""))

# Write output definitions file
write_output_definitions(outputnames)

# save version information
bwa.version <- system(paste(bwa.binary," 2>&1 | grep Version"),intern = TRUE)
documentVersion("BWA",bwa.version)

samtools.version <- system(paste(samtools.binary,"--version | grep samtools"),intern = TRUE)
documentVersion("Samtools",samtools.version)
