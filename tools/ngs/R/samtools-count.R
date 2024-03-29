# TOOL samtools-count.R: "Count alignments in BAM" (Counts alignments in a BAM file, taking the mapping quality into account if needed. This tool is based on the SAMtools package.)
# INPUT alignment.bam TYPE GENERIC
# OUTPUT alignment-counts.txt
# PARAMETER OPTIONAL mapping.quality: "Minimum mapping quality required for counting" TYPE INTEGER FROM 0 TO 1000 DEFAULT 0 (Counts only alignments which have a mapping quality higher than this.)


# EK 26.10.2011

# samtools binary
samtools.binary <- c(file.path(chipster.tools.path, "samtools-1.2", "samtools"))

system(paste(samtools.binary, "view -bc -q", mapping.quality, "alignment.bam > alignment-counts.txt"))
