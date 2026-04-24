# TOOL samtools-merge.R: "Merge BAM" (Merges sorted BAM files and creates an index for the result file. This tool is based on the SAMtools package.)
# INPUT alignment{...}.bam: alignment{...}.bam TYPE GENERIC
# OUTPUT merged.bam
# OUTPUT merged.bam.bai
# RUNTIME R-4.5.1-samtools
# TOOLS_BIN ""

# EK 27.10.2011

# samtools binary
samtools.binary <- c(file.path("/opt/chipster/tools", "samtools", "bin", "samtools"))

# convert sam to bam
system(paste(samtools.binary, "merge merged-not-sorted.bam alignment*.bam"))

# sort bam
system(paste(samtools.binary, "sort merged-not-sorted.bam -o merged.bam"))

# index bam
system(paste(samtools.binary, "index merged.bam"))
