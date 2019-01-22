# TOOL gatk4-genomicsdbiimport.R: "Consolidate GVCFs using GenomicsDBImport with GATK4" (Import single-sample GVCFs into GenomicsDB before joint genotyping. Input is one or more GVCFs produced by in HaplotypeCaller with output type set to GVCF or BP_RESOLUTION, containing the samples to joint-genotype. Currently the tool only supports diploid data. This tool is based on the GATK tool GenomicsDBImport.)
# INPUT variants{...}.g.vcf: "gVCF files" TYPE GENERIC
# OUTPUT OPTIONAL genomicsDBI.tar
# OUTPUT OPTIONAL gatk_log.txt
# PARAMETER OPTIONAL gatk.interval: "Genomic interval" TYPE STRING (A single genomic interval over which to operate. Format chromosome:begin-end, e.g. 20:10,000,000-10,200,000)

# AMS 2019-01-22

source(file.path(chipster.common.path, "tool-utils.R"))
source(file.path(chipster.common.path, "gatk-utils.R"))

# binaries
gatk.binary <- c(file.path(chipster.tools.path, "GATK4", "gatk"))

# Pre-process input files
#
# gVCF file(s)
inputs <- paste("")
input.names <- read.table("chipster-inputs.tsv", header=F, sep="\t")
for (i in 1:nrow(input.names)) {
	if (grepl(".vcf", input.names[i,1])){
		# Index VCF
		formatGatkVcf(input.names[i,1])
		# Add input to command
		inputs <- paste(inputs, "-V", paste(input.names[i,1], ".gz", sep=""))
	}	
}

# Command
command <- paste(gatk.binary, "GenomicsDBImport", inputs)

# Options
if (nchar(gatk.interval) > 0 ){
	command <- paste(command, "-L", gatk.interval)
}

command <- paste(command, "--genomicsdb-workspace-path gdbi")

# Capture stderr
command <- paste(command, "2>> error.txt")

# Run command
system(command)

# Tar ouput folder
system("tar cf genomicsDBI.tar gdbi")

# Return error message if no result
if (!dir.exists("gdbi")){
	system("mv error.txt gatk_log.txt")
}
