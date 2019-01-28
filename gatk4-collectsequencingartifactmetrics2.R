# TOOL gatk4-collectsequencingartifactmetrics.R: "Generate sequencing artifact metrics for filtering by orientation bias with GATK4" (Generate metrics for filtering by orientation bias with GATK4. Tool is based on GATK4  CollectSequencingArtifactMetrics.)
# INPUT tumor.bam: "Tumor BAM" TYPE BAM
# INPUT reference: "Reference genome FASTA" TYPE GENERIC
# OUTPUT OPTIONAL artifact_metrics.txt
# OUTPUT OPTIONAL gatk_log.txt

source(file.path(chipster.common.path, "gatk-utils.R"))
source(file.path(chipster.common.path, "tool-utils.R"))
source(file.path(chipster.common.path, "zip-utils.R"))

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

# Pre-process BAM
system(paste(samtools.binary, "index tumor.bam > tumor.bam.bai"))


# Run CollectSequencingArtifactMetrics
command <- paste(gatk.binary, "CollectSequencingArtifactMetrics", "-O artifact_metrics.txt", "-I tumor.bam", "-R reference.fasta")


runExternal(command)

# Return error message if no result
if (fileNotOk("artifact_metrics.txt")){
	system("mv stderr.log gatk_log.txt")
}

# read input names
inputnames <- read_input_definitions()

basename <- strip_name(inputnames$tumor.bam)

# Make a matrix of output names
outputnames <- matrix(NA, nrow=1, ncol=2)
outputnames[1,] <- c("artifact_metrics.txt", paste(basename, "_artifact_metrics.txt", sep=""))


# Write output definitions file
write_output_definitions(outputnames)
