# TOOL samtools-create-sam-preview.R: "Create a preview for BAM" (Creates a SAM format preview from a BAM file. The preview contains the header and 200 first alignments. This tool is based on the SAMtools package.)
# INPUT alignment.bam TYPE GENERIC
# OUTPUT preview.sam

# samtools binary
samtools.binary <- c(file.path(chipster.tools.path,"samtools-0.1.19","samtools"))

# Print header to preview
system(paste(samtools.binary,"view -H alignment.bam > preview.sam"))
# Print 200 first reads/alignments to preview
system(paste(samtools.binary,"view alignment.bam | head -200 >> preview.sam"))

# Handle output names
source(file.path(chipster.common.path,"tool-utils.R"))

# read input names
inputnames <- read_input_definitions()
base <- strip_name(inputnames$alignment.bam)

# Make a matrix of output names
outputnames <- matrix(NA,nrow = 1,ncol = 2)
outputnames[1,] <- c("preview.sam",paste(base,".prev.sam",sep = ""))

# Write output definitions file
write_output_definitions(outputnames)






