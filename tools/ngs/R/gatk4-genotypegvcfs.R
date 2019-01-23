# TOOL gatk4-genotypegvcfs.R: "Perform joint genotyping on gVCF files with GATK4" (This tool is based on the GATK tool GenotypeGVCFs.)
# INPUT gdb.tar: "GenomicsDB" TYPE GENERIC
# INPUT reference.fasta: "Reference genome" TYPE GENERIC
# OUTPUT OPTIONAL jointcalls.vcf
# OUTPUT OPTIONAL gatk_log.txt
# PARAMETER gatk.ploidy: "Ploidy" TYPE INTEGER DEFAULT 2 (Ploidy per sample. For pooled data, set to (number of samples in each pool * sample ploidy\).)
# PARAMETER OPTIONAL gatk.interval: "Genomic intervals" TYPE STRING (One or more genomic intervals over which to operate. Format chromosome:begin-end, e.g. 20:10,000,000-10,200,000)
# PARAMETER OPTIONAL gatk.padding: "Interval padding" TYPE INTEGER DEFAULT 0 (Amount of padding in bp to add to each interval.)

# AMS 12.07.2016

source(file.path(chipster.common.path, "gatk-utils.R"))
source(file.path(chipster.common.path, "zip-utils.R"))

# binaries
gatk.binary <- c(file.path(chipster.tools.path, "GATK4", "gatk"))

# Pre-process input files
#
# Reference fasta
formatGatkFasta("reference.fasta")
system("mv reference.fasta.dict reference.dict")

# GenomicDB tar file
system("mkdir gdb")
system("tar xf gdb.tar -C gdb --strip-components=1")

# Command
command <- paste(gatk.binary, "GenotypeGVCFs", "-R reference.fasta", "-V gendb://gdb")

# Options
command <- paste(command, "-ploidy", gatk.ploidy)
if (nchar(gatk.interval) > 0 ){
	command <- paste(command, "-L", gatk.interval)
	command <- paste(command, "-ip", gatk.padding)
}

command <- paste(command, "-O jointcalls.vcf")

# Capture stderr
command <- paste(command, "2>> error.txt")

# Run command
system(command)

# Return error message if no result
if (file.exists("jointcalls.vcf") == FALSE){
	system("ls -l >> error.txt")
	system("mv error.txt gatk_log.txt")
}
