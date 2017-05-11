# TOOL single-cell-tag-read-with-gene-exon.R: "Tag read with gene exon" (Adds a BAM tag GE onto reads when the read overlaps the exon of a gene. This tag contains the name of the gene, as reported in the GTF annotations file. )
# INPUT merged.bam: "Merged BAM" TYPE GENERIC
# OUTPUT OPTIONAL merged_tagged.bam
# OUTPUT OPTIONAL log.txt
# PARAMETER OPTIONAL organism: "GTF" TYPE [Homo_sapiens.GRCh38.87, Mus_musculus.GRCm38.87] DEFAULT Homo_sapiens.GRCh38.87 (GTF file to be used in tagging.)


# OUTPUT OPTIONAL log.txt

# ML 12.10.2016 created

path.dropseq <- c(file.path(chipster.tools.path, "drop-seq_tools"))
gtf.path <- "/opt/chipster/genomes/gtf/"

# command start
command <- paste(path.dropseq, "/TagReadWithGeneExon I=merged.bam O=merged_tagged.bam ANNOTATIONS_FILE=", gtf.path, organism, ".gtf TAG=GE  2>> log.txt", sep="")

# run the tool
system(command)

# stop(paste('CHIPSTER-NOTE: ', command))

#EOF