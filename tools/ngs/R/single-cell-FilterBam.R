# TOOL single-cell-FilterBam.R: "Filter and trim BAM" (Remove reads where the cell or molecular barcode has low quality bases, trim adapters and polyA tails.)
# INPUT tagged.bam: "Tagged BAM" TYPE GENERIC
# OUTPUT OPTIONAL unaligned_tagged_polyA_filtered.bam
# OUTPUT OPTIONAL adapter_trim_summary.txt
# OUTPUT OPTIONAL polyA_trimming_report.txt
# PARAMETER OPTIONAL sequence: "Adapter sequence" TYPE STRING DEFAULT AAGCAGTGGTATCAACGCAGAGTGAATGGG (Adapter sequence to trim off. As a default, SMART adapter sequence.)
# PARAMETER OPTIONAL mismatches: "Mismatches in adapter" TYPE INTEGER DEFAULT 0 (How many mismatches allowed in the adapter sequence)
# PARAMETER OPTIONAL num_bases: "Number of bases to check in adapter" TYPE INTEGER DEFAULT 5 (How many bases to check of the adapter sequence)
# PARAMETER OPTIONAL mismatches_polyA: "Mismatches in polyA" TYPE INTEGER DEFAULT 0 (How many mismatches allowed in the polyA sequence)
# PARAMETER OPTIONAL num_bases_polyA: "Number of bases to check in polyA" TYPE INTEGER DEFAULT 6 (How many bases to at least have to be A in the polyA tail)


# PARAMETER OPTIONAL tag_reject: "TAG to reject" TYPE STRING DEFAULT XQ (Tag in the reads that were below quality threshold)
# OUTPUT OPTIONAL log.txt


# ML 08.11.2016 created

path.dropseq <- c(file.path(chipster.tools.path, "drop-seq_tools"))

# FilterBAM:
command <- paste(path.dropseq, "/FilterBAM TAG_REJECT=XQ INPUT=tagged.bam OUTPUT=unaligned_tagged_filtered.bam  2>> log.txt", sep="")
system(command)


# TrimStartingSequence:
# command start
command.start <- paste(path.dropseq, "/TrimStartingSequence INPUT=unaligned_tagged_filtered.bam OUTPUT=unaligned_tagged_trimmed.bam OUTPUT_SUMMARY=adapter_trim_summary.txt", sep="")
# parameters
command.parameters <- paste("SEQUENCE=", sequence, "MISMATCHES=", mismatches, "NUM_BASES=", num_bases)
# run the tool
command <- paste(command.start, command.parameters, " 2>> log.txt")
system(command)


# PolyATrimmer:
# command start
command.start <- paste(path.dropseq, "/PolyATrimmer INPUT=unaligned_tagged_trimmed.bam OUTPUT=unaligned_tagged_polyA_filtered.bam OUTPUT_SUMMARY=polyA_trimming_report.txt", sep="")
# parameters
command.parameters <- paste("MISMATCHES=", mismatches_polyA, "NUM_BASES=", num_bases_polyA)
# run the tool
command <- paste(command.start, command.parameters, " 2>> log.txt")
system(command)


# stop(paste('CHIPSTER-NOTE: ', command))

#EOF
