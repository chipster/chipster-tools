# TOOL gatk4-getpileupsummaries.R: "GATK4 -Tabulate pileup metrics for inferring contamination" (Summarizes counts of reads that support reference, alternate and other alleles for given sites. Results can be used with CalculateContamination. Tool is based on GATK4 GetPileupSummaries.)
# INPUT reads.bam: "BAM file" TYPE BAM
# INPUT variants.vcf: "Variants VCF" TYPE GENERIC
# INPUT OPTIONAL intervals.vcf: "Intervals VCF" TYPE GENERIC
# INPUT OPTIONAL reference: "Reference genome" TYPE GENERIC
# OUTPUT OPTIONAL GetPileupSummaries.tsv
# OUTPUT OPTIONAL gatk_log.txt
# PARAMETER organism: "Reference sequence" TYPE [other, "FILES genomes/fasta .fa"] DEFAULT other (Reference sequence.)
# PARAMETER chr: "Chromosome names in my BAM file look like" TYPE [chr1, 1] DEFAULT 1 (Chromosome names must match in the BAM file and in the reference sequence. Check your BAM and choose accordingly. This only applies to provided reference genomes.)
# PARAMETER OPTIONAL gatk.interval: "Genomic intervals" TYPE STRING (One or more genomic intervals over which to operate. Format chromosome:begin-end, e.g. 20:10,000,000-10,200,000)
# PARAMETER OPTIONAL gatk.padding: "Interval padding" TYPE INTEGER DEFAULT 0 (Amount of padding in bp to add to each interval.)
# PARAMETER OPTIONAL usevariants: "Use Variants VCF also for intervals" TYPE [yes, no] DEFAULT no (Do you wish to use same VCF file for both variants and intervals.)

source(file.path(chipster.common.path, "gatk-utils.R"))
source(file.path(chipster.common.path, "tool-utils.R"))
source(file.path(chipster.common.path, "zip-utils.R"))
source(file.path(chipster.common.path, "vcf-utils.R"))

unzipIfGZipFile("reference")

# binaries
gatk.binary <- c(file.path(chipster.tools.path, "GATK4", "gatk"))
samtools.binary <- c(file.path(chipster.tools.path, "samtools-0.1.19", "samtools"))

# If user provided fasta we use it, else use internal fasta
if (organism == "other"){
	# If user has provided a FASTA, we use it
	if (file.exists("reference")){
		file.rename("reference", "reference.fasta")
	}
}else{
	# If not, we use the internal one.
	internal.fa <- file.path(chipster.tools.path, "genomes", "fasta", paste(organism,".fa",sep="",collapse=""))
	# If chromosome names in BAM have chr, we make a temporary copy of fasta with chr names, otherwise we use it as is.
	if(chr == "chr1"){
		source(file.path(chipster.common.path, "seq-utils.R"))
		addChrToFasta(internal.fa, "reference.fasta") 
	}else{
		file.copy(internal.fa, "reference.fasta")
	}
}	

# Pre-process input files
#
# FASTA
options <- ""
if (fileOk("reference.fasta")){
	formatGatkFasta("reference.fasta")
	system("mv reference.fasta.dict reference.dict")	
	options <- paste(options, "--reference reference.fasta")
}
# BAM
system(paste(samtools.binary, "index reads.bam > reads.bam.bai"))
options <- paste(options, "-I reads.bam")
# VCF
if (fileOk("variants.vcf")){
	formatGatkVcf("variants.vcf",chr)
	options <-paste(options, "-V variants.vcf.gz")
}
if (fileOk("intervals.vcf")){
	formatGatkVcf("intervals.vcf",chr)
	options <-paste(options, "-L intervals.vcf.gz")	
}else if (usevariants == "yes"){
	options <-paste(options, "-L variants.vcf.gz")
}

if (nchar(gatk.interval) > 0 ){
	command <- paste(command, "-L", gatk.interval)
	if (gatk.padding > 0){
		command <- paste(command, "-ip", gatk.padding)
	}
}


# Run GetPileupSummaries
command <- paste(gatk.binary, "GetPileupSummaries", options, "-O GetPileupSummaries.tsv")

runExternal(command)

# Return error message if no result
if (fileNotOk("GetPileupSummaries.tsv")){
	system("mv stderr.log gatk_log.txt")
}

# read input names
inputnames <- read_input_definitions()

# Make a matrix of output names
outputnames <- matrix(NA, nrow=1, ncol=2)
outputnames[1,] <- c("GetPileupSummaries.tsv", paste(strip_name(inputnames$reads.bam), "_getpileupsummaries.tsv", sep=""))

# Write output definitions file
write_output_definitions(outputnames)
