# TOOL gatk_haplotypecaller.R: "Call germline SNPs and indels via local re-assembly of haplotypes with GATK" (HaplotypeCaller is capable of calling SNPs and indels simultaneously via local de-novo assembly of haplotypes in an active region. This tool is based on the GATK package. Please note GATK is licensed for acadeioc use only.)
# INPUT reference.fasta: "Reference genome" TYPE GENERIC
# INPUT alignment{...}.bam: "BAM files" TYPE BAM
# OUTPUT OPTIONAL variants.vcf
# OUTPUT OPTIONAL variants.g.vcf
# OUTPUT OPTIONAL haplotypes.bam
# OUTPUT OPTIONAL gatk_log.txt
# PARAMETER gatk.ploidy: "Ploidy" TYPE INTEGER DEFAULT 2 (Ploidy (number of chromosomes\) per sample. For pooled data, set to (Number of samples in each pool * Sample Ploidy\).)
# PARAMETER OPTIONAL gatk.interval: "genomic intervals" TYPE STRING (One or more genomic intervals over which to operate. Format chromosome:begin-end, e.g. 20:10,000,000-10,200,00)
# PARAMETER OPTIONAL gatk.padding: "Interval padding" TYPE INTEGER DEFAULT 0 (Amount of padding in bp to add to each interval.)
# PARAMETER OPTIONAL gatk.pcrmodel: "PCR indel model to use" TYPE [NONE, HOSTILE, AGGRESSIVE, CONSERVATIVE] DEFAULT NONE (NONE: no specialized PCR error model will be applied; if base insertion/deletion qualities are present they will be used, HOSTILE: a most aggressive model will be applied that sacrifices true positives in order to remove more false positives, AGGRESSIVE: a more aggressive model will be applied that sacrifices true positives in order to remove more false positives, CONSERVATIVE: a less aggressive model will be applied that tries to maintain a high true positive rate at the expense of allowing more false positives.)
# PARAMETER OPTIONAL gatk.bamout: "Output assembled haplotypes" TYPE [yes, no] DEFAULT no (The assembled haplotypes and locally realigned reads will be written as BAM. This is intended to be used only for troubleshooting purposes, in specific areas where you want to better understand why the caller is making specific calls. Turning on this mode may result in serious performance cost for the caller.)
# PARAMETER OPTIONAL gatk.forceactive: "Tag all bases as active" TYPE [yes, no] DEFAULT no (If selected, the active region walker treats all bases as active.)
# PARAMETER OPTIONAL gatk.disableoptimizations: "Disable optimizations" TYPE [yes, no] DEFAULT no (If set, certain early exit optimizations in HaplotypeCaller will be disabled. They aim to save compute and time by skipping calculations if an ActiveRegion is determined to contain no variants.)
# PARAMETER OPTIONAL gatk.erc: "Mode for emitting reference confidence scores" TYPE [NONE, BP_RESOLUTION, GVCF] DEFAULT NONE (Mode for emitting reference confidence scores. NONE: Regular calling without emitting reference confidence calls, BP_RESOLUTION: Reference model emitted site by site, GVCF: Reference model emitted with condensed non-variant blocks, i.e. the GVCF format.)

# AMS 11.07.2016

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
# Create dictionary file
system(paste("java -jar", picard.binary, "CreateSequenceDictionary R=reference.fasta O=reference.dict"))
# BAM file(s)
inputs <- paste("")
input.names <- read.table("chipster-inputs.tsv", header=F, sep="\t")
for (i in 1:nrow(input.names)) {
	if (grepl(".bam", input.names[i,1])){
		# Index BAM
		index.name <- paste(input.names[i,1], ".bai", sep="")
		system(paste(samtools.binary, "index", input.names[i,1], ">", index.name))
		# Add input to command
		inputs <- paste(inputs, "-I", input.names[i,1])
	}	
}

# Command
command <- paste("java -jar", gatk.binary, "-T HaplotypeCaller", "-R reference.fasta", inputs)

# Options
command <- paste(command, "-ploidy", gatk.ploidy)
if (nchar(gatk.interval) > 0 ){
	command <- paste(command, "-L", gatk.interval)
	command <- paste(command, "-ip", gatk.padding)
}
command <- paste(command, "-pcrModel", gatk.pcrmodel)
if (gatk.bamout == "yes"){
	command <- paste(command, "-bamout haplotypes.bam") 
}
if (gatk.forceactive == "yes"){
	command <- paste(command, "-forceActive") 
}
if (gatk.disableoptimizations == "yes"){
	command <- paste(command, "-disableOptimizations") 
}
command <- paste(command, "-ERC", gatk.erc)
# Output file name has to be *.g.vcf for GVCF option to funtion automaticly
if (gatk.erc == "NONE"){
	command <- paste(command, "-o variants.vcf")
}else{
	command <- paste(command, "-o variants.g.vcf")
}

# Capture stderr
command <- paste(command, "2>> error.txt")

# Run command
system(command)

# Return error message if no result
if ((file.exists("variants.vcf") == FALSE) && (file.exists("variants.g.vcf") == FALSE)){
	system("mv error.txt gatk_log.txt")
}
