# TOOL samtools-alignment-statistics.R: "Count alignment statistics for BAM" (Counts statistics for alignments in a BAM file. This tool is based on the SAMtools flagstat tool.)
# INPUT alignment.bam TYPE GENERIC
# OUTPUT alignment-statistics.txt
# RUNTIME R-4.5.1-samtools
# TOOLS_BIN ""

# AMS 17.09.2015

# samtools binary
samtools.binary <- c(file.path("/opt/chipster/tools", "samtools", "bin", "samtools"))

system(paste(samtools.binary, "flagstat", "alignment.bam > alignment-statistics.txt"))
