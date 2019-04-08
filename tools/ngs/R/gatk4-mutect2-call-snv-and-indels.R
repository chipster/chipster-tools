# TOOL gatk4-mutect2-call-snv-and-indels.R: "GATK4 -Call somatic SNVs and INDELs with Mutect2" (Call somatic short variants via local assembly of haplotypes. Short variants include single nucleotide (SNV\) and insertion and deletion (indel\) variants. Tool is based on GATK4 Mutect2 tool.)
# INPUT sample.bam: "Tumor BAM file" TYPE BAM
# INPUT OPTIONAL control.bam: "Normal BAM file" TYPE BAM
# INPUT reference: "Reference genome FASTA" TYPE GENERIC
# INPUT OPTIONAL germline_resource.vcf: "Germline resource VCF" TYPE GENERIC
# INPUT OPTIONAL normal_panel.vcf: "Panel of Normals" TYPE GENERIC
# INPUT OPTIONAL gatk_interval.vcf: "Intervals VCF" TYPE GENERIC
# OUTPUT OPTIONAL mutect2.vcf
# OUTPUT OPTIONAL gatk_log.txt
# OUTPUT OPTIONAL mutect2.bam
# PARAMETER organism: "Reference sequence" TYPE [other, "FILES genomes/fasta .fa"] DEFAULT other (Reference sequence.)
# PARAMETER chr: "Chromosome names in my BAM file look like" TYPE [chr1, 1] DEFAULT 1 (Chromosome names must match in the BAM file and in the reference sequence. Check your BAM and choose accordingly. This only applies to provided reference genomes.)
# PARAMETER tumor: "Sample sample name" TYPE STRING (BAM sample name of sample.)
# PARAMETER normal: "Control sample name" TYPE STRING (BAM sample name of control.)
# PARAMETER OPTIONAL gatk.interval: "Genomic intervals" TYPE STRING (One or more genomic intervals over which to operate. Format chromosome:begin-end, e.g. 20:10,000,000-10,200,000)
# PARAMETER OPTIONAL gatk.padding: "Interval padding" TYPE INTEGER DEFAULT 0 (Amount of padding in bp to add to each interval.)
# PARAMETER OPTIONAL gatk.bamout: "Output assembled haplotypes as BAM" TYPE [yes, no] DEFAULT no (Output assembled haplotypes as BAM.)


## PARAMETER gatk.afofalleles: "Allele fraction of alleles not in germline resource" TYPE DECIMAL DEFAULT -1 (Population allele fraction assigned to alleles not found in germline resource. Only applicable if germline resource file is provided. -1 = use default value. Default for case-only calling is 5e-8 and for matched-control calling 1e-5.)


source(file.path(chipster.common.path, "gatk-utils.R"))
source(file.path(chipster.common.path, "tool-utils.R"))
source(file.path(chipster.common.path, "zip-utils.R"))

unzipIfGZipFile("reference")

# binaries
gatk.binary <- c(file.path(chipster.tools.path, "GATK4", "gatk"))
samtools.binary <- c(file.path(chipster.tools.path, "samtools", "samtools"))

# If user provided fasta we use it, else use internal fasta
if (organism == "other"){
	# If user has provided a FASTA, we use it
	if (file.exists("reference")){
		file.rename("reference", "reference.fasta")
		formatGatkFasta("reference.fasta")
		system("mv reference.fasta.dict reference.dict")
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
# BAM file(s)
system(paste(samtools.binary, "index sample.bam > sample.bam.bai"))
system(paste(samtools.binary, "index control.bam > control.bam.bai"))
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
command <- paste(gatk.binary, "Mutect2", "-O mutect2.vcf", "-R reference.fasta", "-I control.bam", "--normal", normal, "-I sample.bam", "--tumor", tumor)

if (fileOk("germline_resource")){
	command <- paste(command, "--germline_resource germline_resource.vcf.gz")
	#if (gatk.afofalleles != -1) {
	#	command <- paste(command, "--af-of-alleles-not-in-resource", gatk.afofalleles)
	#}
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

if (nchar(gatk.interval) > 0 ){
	command <- paste(command, "-L", gatk.interval)
	if (gatk.padding > 0){
		command <- paste(command, "-ip", gatk.padding)
	}
}

if (gatk.bamout == "yes"){
	command <- paste(command, "-bamout mutect2.bam")
}

# Capture stderr
command <- paste(command, "2>> error.txt")

# Run command
system(command)

# Return error message if no result
if (fileNotOk("mutect2.vcf")){
	system("ls -l >> error.txt")
	system("mv error.txt gatk_log.txt")
}

# read input names
inputnames <- read_input_definitions()

basename <- strip_name(inputnames$sample.bam)

# Make a matrix of output names
outputnames <- matrix(NA, nrow=2, ncol=2)
outputnames[1,] <- c("mutect2.bam", paste(basename, "_mutect2.bam", sep=""))
outputnames[2,] <- c("mutect2.vcf", paste(basename, "_mutect2.vcf", sep=""))

# Write output definitions file
write_output_definitions(outputnames)
