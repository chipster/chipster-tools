# TOOL samtools-alignment-statistics.R: "Count alignment statistics for BAM" (Counts statistics for alignments in a BAM file. This tool is based on the SAMtools flagstat tool.)
# INPUT alignment.bam TYPE GENERIC
# OUTPUT alignment-statistics.txt

# AMS 17.09.2015

# samtools binary
samtools.binary <- c(file.path(chipster.tools.path, "samtools-0.1.19", "samtools"))

system(paste(samtools.binary, "flagstat", "alignment.bam > alignment-statistics.txt"))
