# TOOL single-cell-tag-read-with-gene-exon.R: "Tag read with gene names" (Adds a BAM tag GE onto reads when the read overlaps the exon of a gene. This tag contains the name of the gene, as reported in the GTF annotations file. )
# INPUT merged.bam: "Merged BAM" TYPE GENERIC
# INPUT OPTIONAL own.gtf: "Own GTF file" TYPE GENERIC
# OUTPUT OPTIONAL merged_tagged.bam
# OUTPUT OPTIONAL merged_tagged.bam.bai
# PARAMETER OPTIONAL organism: "GTF" TYPE [other, "FILES genomes/gtf .gtf"] DEFAULT other (GTF file to be used in tagging. No need to select anything here if you are using your own GTF file.)


# OUTPUT OPTIONAL log.txt

# ML 12.10.2016 created
# ML 04.07.2017 added option to use own GTF

## Source required functions
source(file.path(chipster.common.path, "tool-utils.R"))

path.dropseq <- c(file.path(chipster.tools.path, "drop-seq_tools"))


# check out if the file is compressed and if so unzip it
source(file.path(chipster.common.path, "zip-utils.R"))
unzipIfGZipFile("own.gtf")

# command start
# If using own GTF:
if (file.exists("own.gtf")){
	library(tools)
	dir <- getwd()
	dir <- file_path_as_absolute(dir)
	command <- paste(path.dropseq, "/TagReadWithGeneExon I=merged.bam O=merged_tagged.bam ANNOTATIONS_FILE=", dir, "/own.gtf TAG=GE  2>> log.txt", sep="")
	
# if using one of the GTFs available on Chipster:			
}else{
	gtf.path <- "/opt/chipster/genomes/gtf/"
	command <- paste(path.dropseq, "/TagReadWithGeneExon I=merged.bam O=merged_tagged.bam ANNOTATIONS_FILE=", gtf.path, organism, ".gtf TAG=GE  2>> log.txt", sep="")
}

# run the tool
system(command)

# Only index if BAM not empty to prevent returning an empty .bai file
if (fileOk("merged_tagged.bam", minsize=100)){
	# Index BAM
	samtools.binary <- file.path(chipster.tools.path, "samtools", "samtools")
	system(paste(samtools.binary, "index merged_tagged.bam > merged_tagged.bam.bai"))
}
# stop(paste('CHIPSTER-NOTE: ', command))

#EOF