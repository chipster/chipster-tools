# TOOL picard-samtofastq.R: "BAM to FASTQ" (Convert BAM file to FASTQ)
# INPUT input.bam: "BAM file" TYPE GENERIC
# OUTPUT OPTIONAL sam_to_fastq.fastq       


# OUTPUT OPTIONAL log.txt

# 2016-10-31 ML

picard.binary <- file.path(chipster.tools.path, "picard-tools", "picard.jar")

# run
command <- paste("java -Xmx4g -jar", picard.binary, "SamToFastq INPUT=input.bam FASTQ=sam_to_fastq.fastq 2>> log.txt")

#stop(paste('CHIPSTER-NOTE: ', command))

system(command)
