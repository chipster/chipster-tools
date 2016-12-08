# TOOL picard-fastqtosam.R: "FASTQ to BAM" (Converts a FASTQ file to an unaligned BAM or SAM file. This tool extracts read sequences and base qualities from the input FASTQ file and writes them out to a new file in unaligned BAM format.)
# INPUT input.fastq.gz: "FASTQ file" TYPE GENERIC
# INPUT OPTIONAL input2.fastq.gz: "FASTQ file, read2 " TYPE GENERIC
# OUTPUT OPTIONAL fastq_to_bam.bam  


# OUTPUT OPTIONAL log.txt

# 2016-10-31 ML

picard.binary <- file.path(chipster.tools.path, "picard-tools", "picard.jar")

# run
# single end:
if (file.exists("input2.fastq.gz")==FALSE) {
	command <- paste("java -Xmx2g -jar", picard.binary, "FastqToSam F1=input.fastq.gz O=fastq_to_bam.bam SM=something  2>> log.txt")
}

# paired end:
if (file.exists("input2.fastq.gz")==TRUE) {
		command <- paste("java -Xmx2g -jar", picard.binary, "FastqToSam F1=input.fastq.gz F2=input2.fastq.gz O=fastq_to_bam.bam SM=something  2>> log.txt")
	}

#stop(paste('CHIPSTER-NOTE: ', command))

system(command)


