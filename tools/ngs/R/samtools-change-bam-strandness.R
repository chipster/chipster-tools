# TOOL samtools-change-bam-strandness.R: "Reverse XS tags in BAM" (Reverses strandedness tags \(XS\) in BAM. Every XS:A:+ tag is changed to XS:A:- and vice versa. If you have stranded RNA-seq data and suspect that a wrong strandedness parameter value was used when aligning the reads to genome with HISAT2 or TopHat, you can correct the XS tags in the BAM file with this tool. XS tags are used by the Cuffdiff and Cufflinks tools when assembling transcripts. Note that the XS tags are not used by HTseq tools when counting reads per genes.)
# INPUT alignment.bam: "BAM file" TYPE GENERIC
# OUTPUT alignment-edited.bam
# OUTPUT alignment-edited.bam.bai

# AMS 2020-01-21

# samtools binary
samtools.binary <- c(file.path(chipster.tools.path,"samtools-0.1.19","samtools"))

# Change command
system(paste(samtools.binary,"view -h alignment.bam | sed s/XS:A:+/XS:A:P/ | sed s/XS:A:-/XS:A:+/ | sed s/XS:A:P/XS:A:-/ |",samtools.binary,"view -bS - > alignment-edited.bam"))

# index bam
system(paste(samtools.binary,"index alignment-edited.bam"))

# Handle output names
source(file.path(chipster.common.path,"tool-utils.R"))

# read input names
inputnames <- read_input_definitions()
base <- strip_name(inputnames$alignment.bam)

# Make a matrix of output names
outputnames <- matrix(NA,nrow = 2,ncol = 2)
outputnames[1,] <- c("alignment-edited.bam",paste(base,".edited.bam",sep = ""))
outputnames[2,] <- c("alignment-edited.bam.bai",paste(base,".edited.bam.bai",sep = ""))

# Write output definitions file
write_output_definitions(outputnames)