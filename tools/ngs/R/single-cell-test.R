# TOOL single-cell-test.R: "Testing single cell tools" (plaplapla.)
# INPUT unaligned.bam: "Unaligned BAM" TYPE GENERIC
# OUTPUT OPTIONAL unaligned_tagged.bam
# OUTPUT OPTIONAL summary.txt
# OUTPUT OPTIONAL run-summary.txt
# PARAMETER OPTIONAL base_range: "Base range" TYPE STRING DEFAULT 1-12 (Which bases correspond to the cell barcode)
# PARAMETER OPTIONAL base_quality: "Base quality" TYPE INTEGER DEFAULT 10 (Base quality)
# PARAMETER OPTIONAL tag_name: "Tag name" TYPE [XC:XC, XM:XM] DEFAULT XC (Which BAM tag to use for the barcodes)

# ML 12.10.2016 created

## setting up TopHat
#tophat.binary <- c(file.path(chipster.tools.path, "tophat2", "tophat2"))
#path.bowtie <- c(file.path(chipster.tools.path, "bowtie2"))
#path.samtools <- c(file.path(chipster.tools.path, "samtools"))
#set.path <-paste(sep="", "PATH=", path.bowtie, ":", path.samtools, ":$PATH")
#path.bowtie.index <- c(file.path(chipster.tools.path, "genomes", "indexes", "bowtie2", organism))
#path.tophat.index <- c(file.path(chipster.tools.path, "genomes", "indexes", "tophat2", organism))

path.dropseq <- c(file.path(chipster.tools.path, "drop-seq_tools"))

# command start
# command.start <- paste("bash -c '", set.path, tophat.binary)
ommand.start <- paste(path.dropseq, "TagBamWithReadSequenceExtended INPUT=unaligned.bam OUTPUT=unaligned_tagged.bam SUMMARY=summary.txt")

# parameters
#command.parameters <- paste("--bowtie1 -r", mate.inner.distance, "--mate-std-dev", mate.std.dev, "-a", min.anchor.length, "-m", splice.mismatches, "-i", min.intron.length, "-I", max.intron.length, "-g", max.multihits, "--library-type fr-unstranded")
command.parameters <- paste("BASE_RANGE=", base_range, "BARCODED_READ=1, DISCARD_READ=False, TAG_NAME=", tag_name, "NUM_BASES_BELOW_QUALITY=1")

# command ending
#command.end <- paste(path.bowtie.index, reads1, reads2, "2>> tophat.log'")

# run the tool
command <- paste(command.start, command.parameters)

echo.command <- paste("echo '",command ,"' 2>> run.log " )
system(echo.command)
system("echo >> run.log")

system(command)


if (!(file.exists("run.log"))){
	system("mv run.log run-summary.txt")
}


#EOF
