# TOOL bwa-mem.R: "BWA MEM for single or paired end reads" (Aligns reads to genomes using the BWA MEM algorithm. If just one read file is given, then a single end analysis is run. If two read files are given, then mapping is done in paired end mode. If more than two read files are given, you also need to provide a list of filenames of the FASTQ files for each direction. Results are sorted and indexed BAM files, which are ready for viewing in the Chipster genome browser. 
# Note that this BWA tool uses publicly available genomes. If you would like to align reads against your own reference genome, please use the tool \"BWA MEM for single or paired end data with own genome\".)
# INPUT reads{...}.fq: "Reads" TYPE GENERIC
# INPUT OPTIONAL reads1.txt: "List of read 1 files" TYPE GENERIC
# INPUT OPTIONAL reads2.txt: "List of read 2 files" TYPE GENERIC
# OUTPUT OPTIONAL bwa.bam 
# OUTPUT OPTIONAL bwa.bam.bai 
# OUTPUT bwa.log
# PARAMETER organism: "Organism" TYPE ["FILES genomes/indexes/bwa .fa"] DEFAULT "SYMLINK_TARGET genomes/indexes/bwa/default .fa" (Genome that you would like to align your reads against.)
# PARAMETER OPTIONAL minseedlen: "Minimum seed length" TYPE INTEGER DEFAULT 19 (Matches shorter than this will be missed when looking for maximal exact matches or MEMs in the first alignment phase.)
# PARAMETER OPTIONAL bandwith: "Maximum gap length" TYPE INTEGER DEFAULT 100 (Gaps longer than this will not be found. Note also scoring matrix and hit length affect the maximum gap length, in addition to this band width parameter.)
# PARAMETER OPTIONAL matchscore: "Match score" TYPE INTEGER DEFAULT 1 (Score for a matching base.)
# PARAMETER OPTIONAL mismatchscore: "Mismatch penalty" TYPE INTEGER DEFAULT 4 (Penalty for a mismatching base.)
# PARAMETER OPTIONAL gapopen: "Gap opening penalty" TYPE INTEGER DEFAULT 6 (Gap opening penalty.)
# PARAMETER OPTIONAL gapextension: "Gap extension penalty" TYPE INTEGER DEFAULT 1 (Gap extension penalty.)
# PARAMETER OPTIONAL clippenalty: "Penalty for end clipping" TYPE INTEGER DEFAULT 5 (Penalty for 5\'- and 3\'-end clipping. When performing the Smith-Waterman extension of the seed alignments, BWA-MEM keeps track of the best score reaching the end of the read. If this score is larger than the best SW score minus the clipping penalty, clipping will not be applied.) 
# PARAMETER OPTIONAL rgid: "Read group identifier" TYPE STRING (Read group identifier. If you want to add the read group line in the BAM file, you have to give this information.)
# PARAMETER OPTIONAL rgsm: "Sample name for read group" TYPE STRING (The name of the sample sequenced in this read group. Note that you have to fill in also the read group identifier parameter for the read group information to appear in the BAM file.)
# PARAMETER OPTIONAL rgpl: "Platform for read group" TYPE [ none: "Not defined", ILLUMINA, SOLID, LS454, HELICOS, PACBIO] DEFAULT none (Platform\/technology used to produce the read. Note that you have to fill in also the read group identifier parameter for the read group information to appear in the BAM file.)
# PARAMETER OPTIONAL rglb: "Library identifier for read group" TYPE STRING (DNA preparation library identifier. The Mark Duplicates tool uses this field to determine which read groups might contain molecular duplicates, in case the same DNA library was sequenced on multiple lanes. Note that you have to fill in also the read group identifier parameter for the read group information to appear in the BAM file.)

# KM 28.08.2015

source(file.path(chipster.common.path, "tool-utils.R"))
source(file.path(chipster.common.path, "bam-utils.R"))
source(file.path(chipster.common.path, "zip-utils.R"))

# check out if the file is compressed and if so unzip it
input.names <- read.table("chipster-inputs.tsv", header=F, sep="\t")
for (i in 1:nrow(input.names)) {
	unzipIfGZipFile(input.names[i,1])	
}

# bwa binary
bwa.binary <- file.path(chipster.tools.path, "bwa", "bwa mem")
# bwa genome
bwa.genome <- file.path(chipster.tools.path, "genomes", "indexes", "bwa", organism)
# samtools binary
samtools.binary <- c(file.path(chipster.tools.path, "samtools-1.2", "samtools"))

#command.start <- paste("bash -c '", bwa.binary)
command.start <-(bwa.binary)

bwa.parameters <- paste("-M", "-k", minseedlen, "-w", bandwith, "-A", matchscore, "-B", mismatchscore, "-O", gapopen, "-E", gapextension, "-L",  clippenalty )

#Read group definition
if ( nchar(rgid) > 0 ){
	rg.string <-("'@RG")
	rg.string <- paste(rg.string, "\\tID:", rgid, sep="")
	if ( nchar(rgsm) > 0 ){
	   rg.string <- paste(rg.string, "\\tSM:", rgsm, sep="")
    }
	if ( rgpl != "none" ){
    	rg.string <- paste(rg.string, "\\tPL:", rgpl, sep="")
	}
	if ( nchar(rglb) > 0 ){
    	rg.string <- paste(rg.string, "\\tLB:", rglb, sep="")
	}
	rg.string <- paste(rg.string, "'", sep="")
	bwa.parameters <- paste(bwa.parameters,  "-R", rg.string )
}

#Read group error message
if ( nchar(rgid) < 1 ){
	if ( nchar(rgsm) > 0 ){
		stop("CHIPSTER-NOTE: Please define identifier for read group")
	}
	if ( rgpl != "none" ){
		stop("CHIPSTER-NOTE: Please define identifier for read group")
	}
	if ( nchar(rglb) > 0 ){
		stop("CHIPSTER-NOTE: Please define identifier for read group")	
	}
}

# Input files

if (fileOk("reads1.txt")){
	if (fileNotOk("reads2.txt")){
		# Case: One list file -> multiple single-end alignments
		reads1.list <- make_input_list("reads1.txt")
		paired.end <- FALSE
	
	}else{
		# Case: Two list files -> multiple paired-end alignments
		reads1.list <- make_input_list("reads1.txt")
		reads2.list <- make_input_list("reads2.txt")
		paired.end <- TRUE
		if (!(identical(intersect(reads1.list, reads2.list), character(0)))){
			stop(paste('CHIPSTER-NOTE: ', "One or more files is listed in both lists."))	
		}
	}	
}else if (fileNotOk("reads002.fq")){
	# Case: No list file, one fastq -> single single-end alignment
	reads1.list <- paste("reads001.fq")
	paired.end <- FALSE
}else if (fileNotOk("reads003.fq")){
	# Case: No list file, two fastq  -> single paired-end alignment
	in.sorted <- input.names[order(input.names[,2]),]
	reads <- grep("reads", in.sorted[,1], value = TRUE)
	reads1.list <- reads[1]
	reads2.list <- reads[2]
	paired.end <- TRUE
}else{
	# Case: No list files, more than two fastq inputs -> error
	stop(paste('CHIPSTER-NOTE: ', "More than two FASTQ files selected, but no list file(s) provided."))
}

for (i in 1:length(reads1.list)) {
	# command ending
	sam.file <- paste(c(as.character(i), ".sam"), collapse="")
	bam.file <- paste(c(as.character(i), ".bam"), collapse="")
		
	if (paired.end){
		command.end <- paste(bwa.genome, reads1.list[i], reads2.list[i], "1>", sam.file, "2>> bwa.log")
	}else{
		command.end <- paste(bwa.genome, reads1.list[i], "1>", sam.file, "2>> bwa.log")
	}
	# run bwa alignment
	bwa.command <- paste(command.start, bwa.parameters, command.end)
	
	#stop(paste('CHIPSTER-NOTE: ', bwa.command))
	system(bwa.command)
	
	# convert sam to bam
	system(paste(samtools.binary, "view -b", sam.file, "-o", bam.file))
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
for (i in 1:nrow(input.names)) {
	sed.command <- paste("s/", input.names[i,1], "/", input.names[i,2], "/", sep="")
	system(paste("sed -i", sed.command, "bwa.log"))
}


## Handle output names
##
## read input names
#inputnames <- read_input_definitions()
#
## Determine base name
#basename <- strip_name(inputnames$reads.fastq)
#
## Make a matrix of output names
#outputnames <- matrix(NA, nrow=2, ncol=2)
#outputnames[1,] <- c("bwa.bam", paste(basename, ".bam", sep =""))
#outputnames[2,] <- c("bwa.bam.bai", paste(basename, ".bam.bai", sep =""))
#
## Write output definitions file
#write_output_definitions(outputnames)

