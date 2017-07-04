# TOOL single-cell-tag-read-with-gene-exon.R: "Tag read with gene exon" (Adds a BAM tag GE onto reads when the read overlaps the exon of a gene. This tag contains the name of the gene, as reported in the GTF annotations file. )
# INPUT merged.bam: "Merged BAM" TYPE GENERIC
# INPUT OPTIONAL own.gtf: "Own GTF file" TYPE GENERIC
# OUTPUT OPTIONAL merged_tagged.bam
# PARAMETER OPTIONAL organism: "GTF" TYPE [Homo_sapiens.GRCh38.87, Mus_musculus.GRCm38.87] DEFAULT Homo_sapiens.GRCh38.87 (GTF file to be used in tagging. No need to select anything here if you are using your own GTF file.)


# OUTPUT OPTIONAL log.txt

# ML 12.10.2016 created
# ML 04.07.2017 added option to use own GTF

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

# stop(paste('CHIPSTER-NOTE: ', command))

#EOF