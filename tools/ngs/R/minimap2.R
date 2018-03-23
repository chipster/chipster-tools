# TOOL minimap2.R: "Minimap2 for mapping single reads and own genome" (Aligns reads to genomes using the minimap2 algorithm. If just one read file is given, then a single end analysis is run. The input files can be in FASTQ or FASTA format. Results are sorted and indexed BAM files, which are ready for viewing in the Chipster genome browser. 
# Note that this minimap2 tool requires that you have imported the reference genome to Chipster in fasta format.)
# INPUT reads.fq : "Reads" TYPE GENERIC
# INPUT genome.fasta: "Reference genome" TYPE GENERIC
# OUTPUT OPTIONAL minimap2.bam 
# OUTPUT OPTIONAL minimap2.bam.bai 
# OUTPUT minimap2.log 
# PARAMETER OPTIONAL task: "Task type" TYPE [ map-pb: "Map PacBio subreads to a genome", map-ont: "Map Oxford nanopore reads to a genome", splice_2: "Map PacBio Iso-seq or traditional cDNA to reference", splice: "Map Nanopore 2D cDNA-seq data to refrence",  splice_3: "Map Nanopore Direct RNA-seq to refrence",  splice_4: "Mapping against SIRV control reference", asm5: "Aligning assebly to reference genome" ] DEFAULT none (Mapping or aligment task to be performed)
# PARAMETER OPTIONAL rgid: "Read group identifier" TYPE STRING (Read group identifier. If you want to add the read group line in the BAM file, you have to give this information.)
# PARAMETER OPTIONAL rgsm: "Sample name for read group" TYPE STRING (The name of the sample sequenced in this read group. Note that you have to fill in also the read group identifier parameter for the read group information to appear in the BAM file.)
# PARAMETER OPTIONAL rgpl: "Platform for read group" TYPE [ none: "Not defined", ILLUMINA, SOLID, LS454, HELICOS, PACBIO] DEFAULT none (Platform\/technology used to produce the read. Note that you have to fill in also the read group identifier parameter for the read group information to appear in the BAM file.)
# PARAMETER OPTIONAL rglb: "Library identifier for read group" TYPE STRING (DNA preparation library identifier. The Mark Duplicates tool uses this field to determine which read groups might contain molecular duplicates, in case the same DNA library was sequenced on multiple lanes. Note that you have to fill in also the read group identifier parameter for the read group information to appear in the BAM file.)

# KM 6.3.2018


# Handle output names
source(file.path(chipster.common.path, "tool-utils.R"))
# read input names
inputnames <- read_input_definitions()

source(file.path(chipster.common.path, "bam-utils.R"))

# check out if the file is compressed and if so unzip it
# above step is not needed as minimap 2 does this automatically 

# minimap2
#minimap2.binary <- file.path(chipster.tools.path, "minimap2", "minimap2")
minimap2.binary <- file.path("/opt/chipster/tools_local/minimap2-2.9_x64-linux/minimap2")
samtools.binary <- file.path(chipster.tools.path, "samtools-1.2", "samtools")

# User should allways decide the analysis mode
# No default task is selected
task.test <- paste("X", task, "X", sep="")
if ( task.test == "XX" ){
	stop("CHIPSTER-NOTE: Please define the aligment task type.")
}


command.start <-(minimap2.binary)

if ( task == "splice_2"){
	task <- ("splice -uf")
}

if ( task == "splice_3"){
	task <- ("splice -uf -k14")
}

if ( task == "splice_4"){
	task <- ("splice -uf -k14")
}

minimap2.parameters <- paste("-ax", task )

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
	minimap2.parameters <- paste(minimap2.parameters,  "-R", rg.string )
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



command.end <- paste("genome.fasta reads.fq 1> alignment.sam 2>> minimap2.log")


# run minimap2 alignment
minimap2.command <- paste(command.start, minimap2.parameters, command.end)

echo.command <- paste("echo '", minimap2.command, "' > minimap2.log")
system(echo.command)
system(minimap2.command)

#
sam.size <- system("du alignment.sam | cut -f1", intern = TRUE )
if ( sam.size < 10 ){
	system ("ls -l >> minimap2.log")
	echo.command <- paste("echo '",task.test,"'>> minimap2.log")
	system (echo.command)
	system ("echo -------------------------------------- >> minimap2.log")
	system ("echo Alingnment file contains no data or is very small >> minimap2.log")
	system ("echo        >> minimap2.log ")
	system ("echo Content of the resulting sam file: >> minimap2.log")
	system ("cat alignment.sam >> minimap2.log")
}else{
# convert sam to bam
system(paste(samtools.binary, "view -b alignment.sam -o alignment.bam"))

#system("ls -l >> minimap2.log")
# Change file named in BAM header to display names
displayNamesToBAM("alignment.bam")

# sort bam
system(paste(samtools.binary, "sort alignment.bam alignment.sorted"))

# index bam
system(paste(samtools.binary, "index alignment.sorted.bam"))
}

# rename result files
system("mv alignment.sorted.bam minimap2.bam")
system("mv alignment.sorted.bam.bai minimap2.bam.bai")
system("ls -l >> minimap2.log")

# Substitute display names to log for clarity
#displayNamesToFile("minimap2.log")

# Handle output names
#
# read input names
inputnames <- read_input_definitions()

# Determine base name

# Default name if neither special case match. 
# Special cases are for compatibility with older versions of the script
basename <- strip_name(inputnames$reads.fq)


# Make a matrix of output names
outputnames <- matrix(NA, nrow=2, ncol=2)
outputnames[1,] <- c("minimap2.bam", paste(basename, "_minimap2.bam", sep =""))
outputnames[2,] <- c("minimap2.bam.bai", paste(basename, "_minimap2.bam.bai", sep =""))

# Write output definitions file
write_output_definitions(outputnames)


