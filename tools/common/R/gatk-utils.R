# Formats a VCF file to be comatible with GATK. 
# Most GATK tools require VCF files to be bgzip compressed and tabix indexed.
#
# Note that e.g. for input file "example.vcf" you end up with files "example.vcf.gz" and "example.vcf.gz.tbi".
# Remember to use correct filenames in GATK commands.
#
formatGatkVcf <- function(input.vcf){
	source(file.path(chipster.common.path, "zip-utils.R"))
	tabix.binary <- c(file.path(chipster.tools.path, "tabix", "tabix"))
	bgzip.binary <- c(file.path(chipster.tools.path, "tabix", "bgzip"))
	# Uncompress
	unzipIfGZipFile(input.vcf)
	# Bgzip
	system(paste(bgzip.binary, input.vcf))
	# Index
	gz.name <- paste(input.vcf, ".gz", sep="")
	system(paste(tabix.binary, "-p vcf", gz.name))
	
}

# Formast a FASTA format reference sequence to be compatible with GATK
# Most GATK tools require the reference sequence to be indexed and have a dictionary file
#
# Note that e.g for input file "example.fa" you end up with files "example.fa.fai" and "example.fa.dict".
#
formatGatkFasta <- function(input.vcf){
	source(file.path(chipster.common.path, "zip-utils.R"))
	picard.binary <- c(file.path(chipster.tools.path, "picard-tools", "picard.jar"))
	samtools.binary <- c(file.path(chipster.tools.path, "samtools", "samtools"))
	# Uncompress
	unzipIfGZipFile(input.vcf)
	# Index
	system(paste(samtools.binary, "faidx", sequence))
	# Create dictionary file
	system(paste("java -jar", picard.binary, "CreateSequenceDictionary", paste("R=", sequence, sep="", collapse=""), paste("O=", sequence, ".dict", sep="", collapse="")))
}
