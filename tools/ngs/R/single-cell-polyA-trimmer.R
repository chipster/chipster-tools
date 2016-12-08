# TOOL single-cell-polyA-trimmer.R: "PolyA trimmer" (Trim away polyA tails from reads.)
# INPUT unaligned_tagged_trimmed.bam: "Unaligned BAM" TYPE GENERIC
# OUTPUT OPTIONAL unaligned_tagged_polyA_filtered.bam
# OUTPUT OPTIONAL polyA_trimming_report.txt
# OUTPUT OPTIONAL log.txt
# PARAMETER OPTIONAL mismatches: "Mismatches" TYPE INTEGER DEFAULT 0 (How many mismatches allowed)
# PARAMETER OPTIONAL num_bases: "Number of bases" TYPE INTEGER DEFAULT 6 (How many bases to check)

# ML 12.10.2016 created


path.dropseq <- c(file.path(chipster.tools.path, "drop-seq_tools"))

# command start
command.start <- paste(path.dropseq, "/PolyATrimmer INPUT=unaligned_tagged_trimmed.bam OUTPUT=unaligned_tagged_polyA_filtered.bam OUTPUT_SUMMARY=polyA_trimming_report.txt", sep="")

# parameters
command.parameters <- paste("MISMATCHES=", mismatches, "NUM_BASES=", num_bases)

# run the tool
command <- paste(command.start, command.parameters, " 2>> log.txt")

system(command)

# stop(paste('CHIPSTER-NOTE: ', command))

#EOF