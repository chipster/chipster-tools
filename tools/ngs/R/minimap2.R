# TOOL minimap2.R: "Minimap2 for mapping reads to genomes" (This tool aligns reads to genomes using the Minimap2 algorithm. Minimap2 is deverloped especially for long reads. If just one input file is given, then the tool assumes that it contains the reads to be aligned and that the refreces genome is selected with the Genome parameter. If two input files are given the other one will be used as the file contang query sequences and the other as containg the reference sequeces.  The input reads file can be in FASTQ or FASTA format. Reference genome should be in fasta format.. Results are sorted and indexed BAM files. .)
# INPUT reads.fq : "Reads" TYPE GENERIC
# INPUT OPTIONAL genome.fasta: "Reference genome" TYPE GENERIC
# OUTPUT OPTIONAL minimap2.bam 
# OUTPUT OPTIONAL minimap2.bam.bai 
# OUTPUT OPTIONAL minimap2.log 
# PARAMETER OPTIONAL chipster_genome: "Genome" TYPE ["FILES genomes/fasta .fa"] DEFAULT "SYMLINK_TARGET genomes/indexes/bowtie2/default .fa" (Genome that you would like to align your reads against. This parameter is ignored if you provide your own refrence genome as the second input file.)
# PARAMETER OPTIONAL task: "Task type" TYPE [ map-pb: "Map PacBio subreads to a genome", map-ont: "Map Oxford nanopore reads to a genome", splice_2: "Map PacBio Iso-seq or traditional cDNA to reference", splice: "Map Nanopore 2D cDNA-seq data to reference",  splice_3: "Map Nanopore Direct RNA-seq to reference",  splice_4: "Mapping against SIRV control reference", asm5: "Aligning assembly to reference genome" ] DEFAULT none (Mapping or aligment task to be performed)
# PARAMETER OPTIONAL rgid: "Read group identifier" TYPE STRING (Read group identifier. If you want to add the read group line in the BAM file, you have to give this information.)
# PARAMETER OPTIONAL rgsm: "Sample name for read group" TYPE STRING (The name of the sample sequenced in this read group. Note that you have to fill in also the read group identifier parameter for the read group information to appear in the BAM file.)
# PARAMETER OPTIONAL rgpl: "Platform for read group" TYPE [ none: "Not defined", ILLUMINA, SOLID, LS454, HELICOS, PACBIO] DEFAULT none (Platform\/technology used to produce the read. Note that you have to fill in also the read group identifier parameter for the read group information to appear in the BAM file.)
# PARAMETER OPTIONAL rglb: "Library identifier for read group" TYPE STRING (DNA preparation library identifier. The Mark Duplicates tool uses this field to determine which read groups might contain molecular duplicates, in case the same DNA library was sequenced on multiple lanes. Note that you have to fill in also the read group identifier parameter for the read group information to appear in the BAM file.)
# PARAMETER OPTIONAL save_log: "Collect a log file" TYPE [yes: Yes, no: No] DEFAULT yes (Collect a log file about the Mimimap2 mapping process.)

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
#minimap2.binary <- file.path("chipster.tools.path,/opt/chipster/tools_local/minimap2-2.9_x64-linux/minimap2")

conda.path <- file.path( chipster.tools.path,"miniconda3","conda_execute")
conda.env <- ("chipster_tools")
conda.tool <- ("minimap2")
conda.def <- paste(conda.env, "/", conda.tool, sep="")
minimap2.binary <- paste(conda.path, conda.def )


samtools.binary <- file.path(chipster.tools.path, "samtools-1.2", "samtools")

# User should allways decide the analysis mode
# No default task is selected
task.test <- paste("X", task, "X", sep="")
if ( task.test == "XX" ){
	stop("CHIPSTER-NOTE: Please define the aligment task type.")
}


#copy genome file if needed
if(!file.exists("genome.fasta")){
	fasta.genome.path <- c(file.path(chipster.tools.path, "genomes", "fasta"))
	genomecopy.command <- paste("cp ",fasta.genome.path, "/",chipster_genome ,".fa  genome.fasta", sep="" )
	echo.command <- paste("echo Using genome ",fasta.genome.path, "/",chipster_genome ,".fa >> minimap2.log ", sep=""  )
	system(echo.command)
	system(genomecopy.command)
	
}

command.start <-(minimap2.binary)

if ( task == "splice_2"){
	task <- ("splice -uf")
}

if ( task == "splice_3"){
	task <- ("splice -uf -k14")
}

if ( task == "splice_4"){
	task <- ("splice --splice-flank=no")
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

genome.test <- paste("X", task, "X", sep="")
if ( task.test == "XX" ){
	stop("CHIPSTER-NOTE: Please define the aligment task type.")
}

command.end <- paste("genome.fasta reads.fq 1> alignment.sam 2>> minimap2.log")
# run minimap2 alignment
minimap2.command <- paste(command.start, minimap2.parameters, command.end)

echo.command <- paste("echo '", minimap2.command, "' >> minimap2.log")
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


if ( save_log == "no") {
	system ("rm -f minimap2.log")
}


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


