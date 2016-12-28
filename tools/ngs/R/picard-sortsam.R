# TOOL picard-sortsam.R: "Sort BAM" (Sort BAM file based on the quaryname.)
# INPUT aligned.sam: "BAM file" TYPE GENERIC
# OUTPUT OPTIONAL aligned_sorted.bam     


# OUTPUT OPTIONAL log.txt


# 2016-10-31 ML

picard.binary <- file.path(chipster.tools.path, "picard-tools", "picard.jar")

# run
command <- paste("java -Xmx2g -jar", picard.binary, "SortSam I=aligned.sam O=aligned_sorted.bam SO=queryname 2>> log.txt")
system(command)
