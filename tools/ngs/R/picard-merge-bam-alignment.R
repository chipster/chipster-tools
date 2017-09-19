# TOOL picard-merge-bam-alignment.R: "Merge aligned and unaligned BAM" (Merge sorted BAM alignment and unaligned, tagged BAM file. Make sure the input files are assigned correctly!)
# INPUT unmapped.bam: "Unmapped BAM" TYPE GENERIC
# INPUT aligned.bam: "Aligned BAM" TYPE GENERIC
# OUTPUT OPTIONAL merged.bam     
# PARAMETER OPTIONAL reference: "Reference" TYPE ["FILES genomes/fasta .fa"] DEFAULT "SYMLINK_TARGET genomes/fasta/default .fa" (Use same reference as in the alignment!)



# OUTPUT OPTIONAL log.txt

# 2016-10-31 ML
# 2017-05-08 ML add sort BAM to this tool
# 2017-07-19 ML Naming of outputs according to the input names


# Handle output names
# Source read_input_definitions and strip_name functions
source(file.path(chipster.common.path, "tool-utils.R"))
# read input names and strip file extension
inputnames <- read_input_definitions()
input1name <- inputnames$aligned.bam # aligned.bam or unmapped.bam for the name?
input1namestripped <-strip_name(input1name)
#write the input file name into log
write(input1namestripped, file = "log.txt")

# Make a matrix of output names
# These override the default ones
outputnames <- matrix(NA, nrow=1, ncol=2)
outputnames[1,] <- c("merged.bam", paste(input1namestripped, "_merged.bam", sep = ""))
#outputnames[2,] <- c("tagging_and_trimming_summary.txt", "tagging_and_trimming_summary.txt")
#outputnames[3,] <- c("tagging_and_trimming_histograms.pdf", "tagging_and_trimming_histograms.pdf")
#outputnames[4,] <- c("unaligned_tagged.bam", paste(input1namestripped, ".bam", sep =  ""))


# Write output definitions file
write_output_definitions(outputnames)


picard.binary <- file.path(chipster.tools.path, "picard-tools", "picard.jar")
path.reference <- "/opt/chipster/genomes/fasta" #/Homo_sapiens.GRCh38.fa"
#path.reference <- "/opt/chipster/genomes/fasta/Homo_sapiens.GRCh38.fa"
#path.tophat.index <- c(file.path(chipster.tools.path, "genomes", "indexes", "tophat2", organism))

# create symlink

#command <- paste("ln -s ", path.reference, "/Homo_sapiens.GRCh38.fa Homo_sapiens.GRCh38.fa", sep="" )
command <- paste("ln -s ", path.reference, "/", reference, ".fa ", reference, ".fa", sep="" )
system(command)		

# create dictionary (dict) (takes less than minute)
#command <- paste("java -Xmx2g -jar ", picard.binary, " CreateSequenceDictionary R=Homo_sapiens.GRCh38.fa O=Homo_sapiens.GRCh38.dict 2>> log.txt", sep="")
command <- paste("java -Xmx2g -jar ", picard.binary, " CreateSequenceDictionary R=", reference, ".fa O=" ,reference, ".dict 2>> log.txt", sep="")
system(command)


# Sort BAM (Picard):
command <- paste("java -Xmx2g -jar", picard.binary, "SortSam I=aligned.bam O=aligned_sorted.bam SO=queryname 2>> log.txt")
system(command)


# Merge files (Picard):
command <- paste("java -Xmx2g -jar ", picard.binary, " MergeBamAlignment UNMAPPED_BAM=unmapped.bam ALIGNED_BAM=aligned_sorted.bam O=merged.bam R=" ,reference, ".fa INCLUDE_SECONDARY_ALIGNMENTS=false PAIRED_RUN=false 2>> log.txt", sep="") 
system(command)
# stop(paste('CHIPSTER-NOTE: ', command))


# EOF