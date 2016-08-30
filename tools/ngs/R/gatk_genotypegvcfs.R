# TOOL gatk_genotypegvcfs.R: "Perform joint genotyping on gVCF files with GATK" (This tool is based on the GATK package. Please note GATK is licensed for academic use only.)
# INPUT alignment{...}.g.vcf: "gVCF files" TYPE GENERIC
# INPUT reference.fasta: "Reference genome" TYPE GENERIC
# OUTPUT OPTIONAL jointcalls.vcf
# OUTPUT OPTIONAL gatk_log.txt
# PARAMETER gatk.ploidy: "Ploidy" TYPE INTEGER DEFAULT 2 (Ploidy (number of chromosomes\) per sample. For pooled data, set to (Number of samples in each pool * Sample Ploidy\).)
# PARAMETER OPTIONAL gatk.interval: "genomic intervals" TYPE STRING (One or more genomic intervals over which to operate. Format chromosome:begin-end, e.g. 20:10,000,000-10,200,000)
# PARAMETER OPTIONAL gatk.padding: "Interval padding" TYPE INTEGER DEFAULT 0 (Amount of padding in bp to add to each interval.)

# AMS 12.07.2016

source(file.path(chipster.common.path, "zip-utils.R"))
unzipIfGZipFile("reference.fasta")

# binaries
gatk.binary <- c(file.path(chipster.tools.path, "GATK", "GenomeAnalysisTK.jar"))
samtools.binary <- c(file.path(chipster.tools.path, "samtools", "samtools"))
picard.binary <- c(file.path(chipster.tools.path, "picard-tools", "picard.jar"))

# Check if GATK is installed
if (file.exists(gatk.binary) == FALSE){
	source(file.path(chipster.common.path, "gatk-utils.R"))
	message <- noGatkMessage()
	stop(paste('CHIPSTER-NOTE: ', message))
}

# Pre-process input files
#
# Index fasta
system(paste(samtools.binary, "faidx reference.fasta"))
# Create dict file
system(paste("java -jar", picard.binary, "CreateSequenceDictionary R=reference.fasta O=reference.dict"))
# gVCF file(s)
inputs <- paste("")
input.names <- read.table("chipster-inputs.tsv", header=F, sep="\t")
for (i in 1:nrow(input.names)) {
	if (grepl(".vcf", input.names[i,1])){
		# Add input to command
		inputs <- paste(inputs, "-V", input.names[i,1])
	}	
}

# Command
command <- paste("java -jar", gatk.binary, "-T GenotypeGVCFs", "-R reference.fasta", inputs)

# Options
command <- paste(command, "-ploidy", gatk.ploidy)
if (nchar(gatk.interval) > 0 ){
	command <- paste(command, "-L", gatk.interval)
	command <- paste(command, "-ip", gatk.padding)
}

command <- paste(command, "-o jointcalls.vcf")

# Capture stderr
command <- paste(command, "2>> error.txt")

# Run command
system(command)

# Return error message if no result
if (file.exists("jointcalls.vcf") == FALSE){
	system("mv error.txt gatk_log.txt")
}
