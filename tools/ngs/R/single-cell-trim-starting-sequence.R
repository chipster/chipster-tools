# TOOL single-cell-trim-starting-sequence.R: "Trim starting sequence" (plaplapla.)
# INPUT unaligned_tagged_filtered.bam: "Unaligned BAM" TYPE GENERIC
# OUTPUT OPTIONAL unaligned_tagged_trimmed.bam
# OUTPUT OPTIONAL summary.txt
# OUTPUT OPTIONAL log.txt
# PARAMETER OPTIONAL sequence: "Sequence" TYPE STRING DEFAULT AAGCAGTGGTATCAACGCAGAGTGAATGGG (Sequence to trim off. As a default, SMART adapter sequence.)
# PARAMETER OPTIONAL mismatches: "Mismatches" TYPE INTEGER DEFAULT 0 (How many mismatches allowed)
# PARAMETER OPTIONAL num_bases: "Number of bases" TYPE INTEGER DEFAULT 5 (How many bases to check)


# ML 12.10.2016 created

path.dropseq <- c(file.path(chipster.tools.path, "drop-seq_tools"))

# command start
command.start <- paste(path.dropseq, "/TrimStartingSequence INPUT=unaligned_tagged_filtered.bam OUTPUT=unaligned_tagged_trimmed.bam OUTPUT_SUMMARY=summary.txt", sep="")

# parameters
command.parameters <- paste("SEQUENCE=", sequence, "MISMATCHES=", mismatches, "NUM_BASES=", num_bases)

# run the tool
command <- paste(command.start, command.parameters, " 2>> log.txt")

system(command)

# stop(paste('CHIPSTER-NOTE: ', command))

#EOF
