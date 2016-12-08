# TOOL single-cell-detect-bead-synthesis-erros.R: "Detect bead synthesis errors" (Detect and repair barcode synthesis errors.)
# INPUT input.bam: "BAM file" TYPE GENERIC
# OUTPUT OPTIONAL cleaned.bam
# OUTPUT OPTIONAL synthesis_stats.txt
# OUTPUT OPTIONAL synthesis_stats_summary.txt
# PARAMETER OPTIONAL num.barcodes: "Number of barcodes" TYPE INTEGER DEFAULT 2000 (Roughly 2x the number of cells)
# PARAMETER OPTIONAL primer.sequence: "Sequence" TYPE STRING DEFAULT AAGCAGTGGTATCAACGCAGAGTGAATGGG (Sequence to trim off. As a default, SMART adapter sequence.)


# ML 10.12.2016 created

path.dropseq <- c(file.path(chipster.tools.path, "drop-seq_tools"))

# command 
command <- paste(path.dropseq, "/DetectBeadSynthesisErrors I=input.bam O=cleaned.bam OUTPUT_STATS=synthesis_stats.txt SUMMARY=synthesis_stats_summary.txt NUM_BARCODES=",num.barcodes," PRIMER_SEQUENCE=", primer.sequence, " 2>>log.txt", sep="")

# run the tool
system(command)

# stop(paste('CHIPSTER-NOTE: ', command))

#EOF