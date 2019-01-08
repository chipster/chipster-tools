# TOOL gatk4-estimate-contamination.R: "Estimate contamination with GATK4" (Uses GetPileupSummaries, CalculateContamination.)
# INPUT tumor.bam: "Tumor BAM file" TYPE BAM
# INPUT contamination.vcf.gz: "Variants VCF" TYPE GENERIC
# INPUT OPTIONAL reference: "Reference genome" TYPE GENERIC
# OUTPUT OPTIONAL GetPileupSummaries.tsv
# OUTPUT OPTIONAL CalculateContamination.tsv
# OUTPUT OPTIONAL gatk_log.txt
# PARAMETER organism: "Reference sequence" TYPE [other, "FILES genomes/fasta .fa"] DEFAULT other (Reference sequence.)
# PARAMETER chr: "Chromosome names in my BAM file look like" TYPE [chr1, 1] DEFAULT 1 (Chromosome names must match in the BAM file and in the reference sequence. Check your BAM and choose accordingly. This only applies to provided reference genomes.)

source(file.path(chipster.common.path, "tool-utils.R"))
source(file.path(chipster.common.path, "zip-utils.R"))

unzipIfGZipFile("reference")

# binaries
gatk.binary <- c(file.path(chipster.tools.path, "GATK4", "gatk-launch"))
picard.binary <- c(file.path(chipster.tools.path, "picard-tools", "picard.jar"))
samtools.binary <- c(file.path(chipster.tools.path, "samtools", "samtools"))
tabix.binary <- c(file.path(chipster.tools.path, "tabix", "tabix"))

# If user provided fasta we use it, else use internal fasta
if (organism == "other"){
	# If user has provided a FASTA, we use it
	if (file.exists("reference")){
		file.rename("reference", "reference.fasta")
	}else{
		stop(paste('CHIPSTER-NOTE: ', "You need to provide a FASTA file or choose one of the provided reference genomes."))
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
# Index fasta
system(paste(samtools.binary, "faidx reference.fasta"))
# Create dictionary file
system(paste("java -jar", picard.binary, "CreateSequenceDictionary R=reference.fasta O=reference.dict"))
# BAM file(s)
system(paste(samtools.binary, "index tumor.bam > tumor.bam.bai"))
# Germline resource
#system(paste(gatk.binary, "IndexFeatureFile -F germline_resource"))
system(paste(tabix.binary, "contamination.vcf.gz"))

# Run GetPileupSummaries
command <- paste(gatk.binary, "GetPileupSummaries", "-O GetPileupSummaries.tsv", "-R reference.fasta", "-I tumor.bam", "-V contamination.vcf.gz")
runExternal(command)

# Run CalculateContamination
command <- paste(gatk.binary, "CalculateContamination", "-I GetPileupSummaries.tsv", "-O CalculateContamination.tsv")
runExternal(command)

# Return error message if no result
if (fileNotOk("CalculateContamination.tsv")){
	system("mv stderr.log gatk_log.txt")
}

