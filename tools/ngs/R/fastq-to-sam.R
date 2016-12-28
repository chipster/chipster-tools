# TOOL fastq-to-sam.R: "Transform FASTQ file to SAM file" (plaplapla.)
# INPUT unaligned.fastq: "FASTQ file" TYPE GENERIC
# OUTPUT OPTIONAL unaligned.bam
# OUTPUT OPTIONAL run-summary.txt

# ML 12.10.2016 created


## setting up TopHat
#tophat.binary <- c(file.path(chipster.tools.path, "tophat2", "tophat2"))
#path.bowtie <- c(file.path(chipster.tools.path, "bowtie2"))
#path.samtools <- c(file.path(chipster.tools.path, "samtools"))
#set.path <-paste(sep="", "PATH=", path.bowtie, ":", path.samtools, ":$PATH")
#path.bowtie.index <- c(file.path(chipster.tools.path, "genomes", "indexes", "bowtie2", organism))
#path.tophat.index <- c(file.path(chipster.tools.path, "genomes", "indexes", "tophat2", organism))

path.picard <- c(file.path(chipster.tools.path, "picard-tools-2.6.0"))

# command start
# command.start <- paste("bash -c '", set.path, tophat.binary)
command <- paste("java jvm-args -jar picard.jar", path.picard, "INPUT=unaligned.fastq OUTPUT=unaligned.bam")


# run the tool

echo.command <- paste("echo '",command ,"' 2>> run.log " )
print("Indexing the genome...")
system("echo Indexing the genome... >> run.log")
system(echo.command)
system("echo >> run.log")

system(command)


if (!(file.exists("run.log"))){
	system("mv run.log run-summary.txt")
}
# stop("CHIPSTER-NOTE: loppuun meni.")
rtdhyyj

#EOF
