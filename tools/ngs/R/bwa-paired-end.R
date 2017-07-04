# TOOL bwa-paired-end.R: "BWA-backtrack for paired end reads" (BWA-backtrack aligns paired end reads to genomes with BWA ALN algorithm. If more than two read files are given, you also need to provide a list of filenames of the FASTQ files for each direction. Results are sorted and indexed BAM files, which are ready for viewing in the Chipster genome browser. 
# Note that this BWA tool uses publicly available genomes. If you would like to align reads against your own datasets, please use the tool \"BWA for paired-end reads and own genome\".)
# INPUT reads{...}.fq: "Reads" TYPE GENERIC
# INPUT OPTIONAL reads1.txt: "List of read 1 files" TYPE GENERIC
# INPUT OPTIONAL reads2.txt: "List of read 2 files" TYPE GENERIC
# OUTPUT bwa.bam 
# OUTPUT bwa.bam.bai 
# OUTPUT bwa.log
# PARAMETER organism: "Genome" TYPE ["FILES genomes/indexes/bwa .fa"] DEFAULT "SYMLINK_TARGET genomes/indexes/bwa/default .fa" (Genome or transcriptome that you would like to align your reads against.)
# PARAMETER quality.format: "Quality value format used" TYPE [solexa1_3: "Illumina GA v1.3-1.5", sanger: Sanger] DEFAULT sanger (Note that this parameter is taken into account only if you chose to apply the mismatch limit to the seed region. Are the quality values in the Sanger format (ASCII characters equal to the Phred quality plus 33\) or in the Illumina Genome Analyzer Pipeline v1.3 or later format (ASCII characters equal to the Phred quality plus 64\)? Please see the manual for details. Corresponds to the command line parameter -I.)
# PARAMETER OPTIONAL seed.length: "Seed length" TYPE INTEGER DEFAULT 32 (Number of first nucleotides to be used as a seed. If the seed length is longer than query sequences, then seeding will be disabled. Corresponds to the command line parameter -l) 
# PARAMETER OPTIONAL seed.edit:"Maximum of differences in the seed" TYPE INTEGER DEFAULT 2 (Maximum differences in the seed. Corresponds to the command line parameter -k )
# PARAMETER OPTIONAL total.edit: "Maximum edit distance for the whole read" TYPE DECIMAL DEFAULT 0.04 ( Maximum edit distance if the value is more than one. If the value is between 1 and 0 then it defines the fraction of missing alignments given 2% uniform base error rate. In the latter case, the maximum edit distance is automatically chosen for different read lengths. Corresponds to the command line parameter -n.)
# PARAMETER OPTIONAL num.gaps: "Maximum number of gap openings" TYPE INTEGER DEFAULT 1 (Maximum number of gap openings for one read. Corresponds to the command line parameter -o.)
# PARAMETER OPTIONAL num.extensions: "Maximum number of gap extensions" TYPE INTEGER DEFAULT -1 (Maximum number of gap extensions, -1 for disabling long gaps. Corresponds to the command line parameter -e.)
# PARAMETER OPTIONAL gap.opening: "Gap open penalty " TYPE INTEGER DEFAULT 11 (Gap opening penalty. Corresponds to the command line parameter -O.)
# PARAMETER OPTIONAL gap.extension: "Gap extension penalty " TYPE INTEGER DEFAULT 4 (Gap extension penalty. Corresponds to the command line parameter -E.)
# PARAMETER OPTIONAL mismatch.penalty: "Mismatch penalty threshold" TYPE INTEGER DEFAULT 3 (BWA will not search for suboptimal hits with a score lower than defined. Corresponds to the command line parameter -M.)
# PARAMETER OPTIONAL disallow.gaps: "Disallow gaps in region"  TYPE INTEGER DEFAULT 16 (Disallow a long deletion within the given number of bp towards the 3\'-end. Corresponds to the command line parameter -d.)
# PARAMETER OPTIONAL disallow.indel: "Disallow indels within the given number of bp towards the ends"  TYPE INTEGER DEFAULT 5 (Do not put an indel within the defined value of bp towards the ends. Corresponds to the command line parameter -i.)
# PARAMETER OPTIONAL trim.threshold: "Quality trimming threshold" TYPE INTEGER DEFAULT 0 (Quality threshold for read trimming down to 35bp. Corresponds to the command line parameter -q.)
# PARAMETER OPTIONAL barcode.length: "Barcode length"  TYPE INTEGER DEFAULT 0 (Length of barcode starting from the 5 pime-end. The barcode of each read will be trimmed before mapping. Corresponds to the command line parameter -B.)
# PARAMETER OPTIONAL alignment.no: "Maximum number of hits to report for paired reads" TYPE INTEGER DEFAULT 3 (Maximum number of alignments to output in the XA tag for reads paired properly. If a read has more than the given amount of hits, the XA tag will not be written. Corresponds to the command line parameter bwa sampe -n.)
# PARAMETER OPTIONAL max.discordant: "Maximum number of hits to report for discordant pairs" TYPE INTEGER DEFAULT 10 (Maximum number of alignments to output in the XA tag for disconcordant read pairs, excluding singletons. If a read has more than INT hits, the XA tag will not be written. Corresponds to the command line parameter bwa sampe -N.) 
# PARAMETER OPTIONAL max.insert: "Maximum insert size" TYPE INTEGER DEFAULT 500 (Maximum insert size for a read pair to be considered being mapped properly. This option is only used when there are not enough good alignments to infer the distribution of insert sizes. Corresponds to the command line parameter bwa sampe -a.)
# PARAMETER OPTIONAL max.occurrence: "Maximum occurrences for one end" TYPE INTEGER DEFAULT 100000 (Maximum occurrences of a read for pairing. A read with more occurrneces will be treated as a single-end read. Reducing this parameter helps faster pairing. The default value is 100000. For reads shorter than 30bp, applying a smaller value is recommended to get a sensible speed at the cost of pairing accuracy. Corresponds to the command line parameter bwa sampe -o.)

# KM 26.8.2011
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


# bwa
bwa.binary <- file.path(chipster.tools.path, "bwa", "bwa")
# samtools binary
samtools.binary <- c(file.path(chipster.tools.path, "samtools-1.2", "samtools"))
#bwa.indexes <- file.path(chipster.tools.path, "bwa_indexes")
bwa.genome <- file.path(chipster.tools.path, "genomes", "indexes", "bwa", organism)
command.start <- paste("bash -c '", bwa.binary)

###
# common parameters for bwa runs
###
# mode specific parameters
if (total.edit >= 1) {
	total.edit <- round(total.edit)
}
quality.parameter <- ifelse(quality.format == "solexa1_3", "-I", "")
mode.parameters <- paste("aln", "-t", chipster.threads.max, "-o", num.gaps, "-e", num.extensions, "-d", disallow.gaps, "-i" , disallow.indel , "-l" , seed.length , "-k" , seed.edit , "-O" , gap.opening , "-E" , gap.extension , "-q" , trim.threshold, "-B" , barcode.length , "-M" , mismatch.penalty , "-n" , total.edit , quality.parameter)

# Check input files
if (fileOk("reads1.txt")){
	if (fileNotOk("reads2.txt")){
		# Only one list file provided
		stop(paste('CHIPSTER-NOTE: ', "You must provide a list file for both directions."))	
	}else{
		# Two list files provided
		reads1.list <- make_input_list("reads1.txt")
		reads2.list <- make_input_list("reads2.txt")
		if (!(identical(intersect(reads1.list, reads2.list), character(0)))){
			stop(paste('CHIPSTER-NOTE: ', "One or more files is listed in both lists."))	
		}
	}	
}else{
	if (fileNotOk("reads003.fq")){
		in.sorted <- input.names[order(input.names[,2]),]
		reads <- grep("reads", in.sorted[,1], value = TRUE)
		reads1.list <- reads[1]
		reads2.list <- reads[2]
	}else{
		stop(paste('CHIPSTER-NOTE: ', "When providing more than two FASTQ files, you must also provide a list file for both directions."))	
	}
}

# Run BWA for each input
for (i in 1:length(reads1.list)) {
	
	sai1.file <- paste(c(as.character(i), ".1.sai"), collapse="")
	sai2.file <- paste(c(as.character(i), ".2.sai"), collapse="")
	sam.file <- paste(c(as.character(i), ".sam"), collapse="")
	bam.file <- paste(c(as.character(i), ".bam"), collapse="")
		
	# Run first set
	command.end <- paste(bwa.genome, reads1.list[i], "1>", sai1.file, "2>> bwa.log'")
	bwa.command <- paste(command.start, mode.parameters, command.end)
	system(bwa.command)
	
	# Run second set
	command.end <- paste(bwa.genome, reads2.list[i], "1>", sai2.file, "2>> bwa.log'")
	bwa.command <- paste(command.start, mode.parameters, command.end)
	system(bwa.command)
	
	# sai to sam conversion
	sampe.parameters <- paste("sampe -n", alignment.no, "-a", max.insert, "-o" , max.occurrence , "-N" , max.discordant )
	sampe.end <- paste(bwa.genome, sai1.file, sai2.file, reads1.list[i], reads2.list[i], "1>", sam.file, "2>>bwa.log'" )
	sampe.command <- paste( command.start, sampe.parameters , sampe.end )
	system(sampe.command)
			
	# convert sam to bam
	system(paste(samtools.binary, "view -bS", sam.file, "-o", bam.file))
}

# Join bam files
if (fileOk("2.bam")){
	# more than one bam exists, so join them
	system("ls *.bam > bam.list")
	system(paste(samtools.binary, "merge -b bam.list alignment.bam"))
}else{
	# only one bam, so just rename it
	system("mv 1.bam alignment.bam")
}

# Change file named in BAM header to display names
displayNamesToBAM("alignment.bam")

# sort bam
system(paste(samtools.binary, "sort alignment.bam alignment.sorted"))

# index bam
system(paste(samtools.binary, "index alignment.sorted.bam"))

# rename result files
system("mv alignment.sorted.bam bwa.bam")
system("mv alignment.sorted.bam.bai bwa.bam.bai")

# Substitute display names to log for clarity
displayNamesToFile("bwa.log")

# Handle output names
#
source(file.path(chipster.common.path, "tool-utils.R"))

# read input names
inputnames <- read_input_definitions()

if (fileNotOk("reads003.fq")){
	# Determine base name
	base1 <- strip_name(inputnames$reads001.fq)
	base2 <- strip_name(inputnames$reads002.fq)
	basename <- paired_name(base1, base2)
}else{
	basename <- "bwa_multi"
}
# Make a matrix of output names
outputnames <- matrix(NA, nrow=2, ncol=2)
outputnames[1,] <- c("bwa.bam", paste(basename, ".bam", sep =""))
outputnames[2,] <- c("bwa.bam.bai", paste(basename, ".bam.bai", sep =""))

# Write output definitions file
write_output_definitions(outputnames)