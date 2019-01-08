# TOOL gatk4-mutect2-call-snv-and-indels.R: "Call somatic SNVs and INDELs with Mutect2" (Mutec2, GATK4)
# INPUT reads.bam: "Reads BAM file" TYPE BAM
# INPUT tumor.bam: "Tumor BAM file" TYPE BAM
# INPUT reference: "Reference genome" TYPE GENERIC
# INPUT OPTIONAL germline_resource.vcf: "Germline resource" TYPE GENERIC
# INPUT OPTIONAL normal_panel.vcf: "Normal panel" TYPE GENERIC
# INPUT OPTIONAL gatk_interval.vcf: "VCF file to use for intervals" TYPE GENERIC
# OUTPUT OPTIONAL mutect2.vcf
# OUTPUT OPTIONAL gatk_log.txt
# OUTPUT OPTIONAL mutect2.bam
# PARAMETER organism: "Reference sequence" TYPE [other, "FILES genomes/fasta .fa"] DEFAULT other (Reference sequence.)
# PARAMETER chr: "Chromosome names in my BAM file look like" TYPE [chr1, 1] DEFAULT 1 (Chromosome names must match in the BAM file and in the reference sequence. Check your BAM and choose accordingly. This only applies to provided reference genomes.)
# PARAMETER tumor: "Tumor sample name" TYPE STRING (BAM sample name of tumor.)
# PARAMETER normal: "Normal sample name" TYPE STRING (BAM sample name of normal.)
# PARAMETER OPTIONAL gatk.interval: "Genomic intervals" TYPE STRING (One or more genomic intervals over which to operate. Format chromosome:begin-end, e.g. 20:10,000,000-10,200,000)
# PARAMETER OPTIONAL gatk.padding: "Interval padding" TYPE INTEGER DEFAULT 0 (Amount of padding in bp to add to each interval.)
# PARAMETER OPTIONAL gatk.disableoptimizations: "Disable optimizations" TYPE [true, false] DEFAULT false (If set, don't skip calculations in ActiveRegions with no variants.)
# PARAMETER OPTIONAL gatk.donttrimactiveregions: "Don't trim active regions" TYPE [true, false] DEFAULT false (If specified, we will not trim down the active region from the full region \(active + extension\) to just the active interval for genotyping.)
# PARAMETER OPTIONAL gatk.activeprobabilitythreshold: "Minimum probability for a locus to be considered active" TYPE DECIMAL DEFAULT 0.002 (Minimum probability for a locus to be considered active.)
# PARAMETER OPTIONAL gatk.bamout: "Output assembled haplotypes as BAM" TYPE [yes, no] DEFAULT no (Output assembled haplotypes as BAM.)

source(file.path(chipster.common.path, "gatk-utils.R"))
source(file.path(chipster.common.path, "tool-utils.R"))
source(file.path(chipster.common.path, "zip-utils.R"))

unzipIfGZipFile("reference")

# binaries
gatk.binary <- c(file.path(chipster.tools.path, "GATK4", "gatk-launch"))
picard.binary <- c(file.path(chipster.tools.path, "picard-tools", "picard.jar"))
samtools.binary <- c(file.path(chipster.tools.path, "samtools", "samtools"))
tabix.binary <- c(file.path(chipster.tools.path, "tabix", "tabix"))
bgzip.binary <- c(file.path(chipster.tools.path, "tabix", "bgzip"))

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
system(paste(samtools.binary, "index reads.bam > reads.bam.bai"))
system(paste(samtools.binary, "index tumor.bam > tumor.bam.bai"))
# VCF files. These need to bgzip compressed and tabix indexed
if (fileOk("germline_resource.vcf")){
	formatGatkVcf("germline_resource.vcf")	
}
if (fileOk("normal_panel.vcf")){
	formatGatkVcf("normal_panel.vcf")
}
if (fileOk("gatk_interval.vcf")){
	formatGatkVcf("gatk_interval.vcf")
}

# Command
command <- paste(gatk.binary, "Mutect2", "-threads", chipster.threads.max, "-O mutect2.vcf", "-R reference.fasta", "-I reads.bam", "--normal", normal, "-I tumor.bam", "--tumor", tumor)

if (fileOk("germline_resource")){
	command <- paste(command, "--germline_resource germline_resource.vcf.gz")
}

if (fileOk("normal_panel")){
	command <- paste(command, "--normal_panel normal_panel.vcf.gz")
}

if (nchar(gatk.interval) > 0 ){
	command <- paste(command, "-L", gatk.interval)
}

if (fileOk("gatk_interval.file")){
	command <- paste(command, "-L", "gatk_interval.vcf.gz")
}

if (gatk.padding > 0){
	command <- paste(command, "-ip", gatk.padding)
}

if (gatk.disableoptimizations == "true"){
	command <- paste(command, "--disableOptimizations")
}

if (gatk.donttrimactiveregions == "true"){
	command <- paste(command, "--dontTrimActiveRegions")
}

command <- paste(command, "--activeProbabilityThreshold", gatk.activeprobabilitythreshold)
		
if (gatk.bamout == "yes"){
	command <- paste(command, "-bamout mutect2.bam")
}

# Capture stderr
command <- paste(command, "2>> error.txt")

# Run command
system(command)

# Return error message if no result
system("mv error.txt gatk_log.txt")

