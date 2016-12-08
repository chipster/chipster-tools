# TOOL single-cell-TagBam.R: "Tag BAM" (Tag the reads based on the cell barcodes using XC tag and molecular barcodes using XM tag. The rest of the read is discarded. XQ tag is added to each read to represent the number of bases that have quality scores below the base quality threshold.)
# INPUT unaligned.bam: "Unaligned BAM" TYPE GENERIC
# OUTPUT OPTIONAL unaligned_tagged.bam
# OUTPUT OPTIONAL summary_cell.txt
# OUTPUT OPTIONAL summary_molecular.txt
# PARAMETER base_range_cell: "Base range for cell barcode" TYPE STRING DEFAULT 1-12 (Which bases correspond to the cell barcode)
# PARAMETER base_range_mol: "Base range for molecule barcode" TYPE STRING DEFAULT 13-20 (Which bases correspond to the molecule barcode)
# PARAMETER OPTIONAL base_quality: "Base quality" TYPE INTEGER DEFAULT 10 (Mark in the XQ tag how many bases fall below this threshold.)


# PARAMETER OPTIONAL discard_read: "Discard read" TYPE [True, False] DEFAULT False (Discard the read)
# PARAMETER OPTIONAL tag_name_cell: "Tag name for cell" TYPE [XC:XC, XM:XM] DEFAULT XC (Which BAM tag to use for the barcodes)
# OUTPUT OPTIONAL log.txt


# ML 12.10.2016 created

path.dropseq <- c(file.path(chipster.tools.path, "drop-seq_tools"))

# First round: cell barcode

# command start
command.start <- paste(path.dropseq, "/TagBamWithReadSequenceExtended INPUT=unaligned.bam OUTPUT=unaligned_tagged_cell.bam SUMMARY=summary_cell.txt", sep="")
# parameters
command.parameters <- paste("BASE_RANGE=", base_range_cell, " BASE_QUALITY=",base_quality ," BARCODED_READ=1 DISCARD_READ=False TAG_NAME=XC NUM_BASES_BELOW_QUALITY=1")
# run the tool
command <- paste(command.start, command.parameters, " 2>> log.txt")
system(command)

# Second round: molecule barcode

# command start
command.start <- paste(path.dropseq, "/TagBamWithReadSequenceExtended INPUT=unaligned_tagged_cell.bam OUTPUT=unaligned_tagged.bam SUMMARY=summary_molecular.txt", sep="")
# parameters
command.parameters <- paste("BASE_RANGE=", base_range_mol, "BASE_QUALITY=", base_quality ," BARCODED_READ=1 DISCARD_READ=True TAG_NAME=XM NUM_BASES_BELOW_QUALITY=1")

# run the tool
command <- paste(command.start, command.parameters, " 2>> log.txt")
system(command)

# stop(paste('CHIPSTER-NOTE: ', command))

#EOF
