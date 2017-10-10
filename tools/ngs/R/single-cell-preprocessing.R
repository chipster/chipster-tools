# TOOL single-cell-preprocessing.R: "Preprocessing DropSeq FASTQ files" (Transforms FASTQ to BAM, then tags the reads based on the cell barcodes using XC tag and molecular barcodes using XM tag. The rest of the read is discarded. XQ tag is added to each read to represent the number of bases that have quality scores below the base quality threshold. Filtering can be performed. Finally the BAM is transformed back to FASTQ format.)
# INPUT input.fastq.gz: "FASTQ file, read1, barcode read" TYPE GENERIC
# INPUT OPTIONAL input2.fastq.gz: "FASTQ file, read2, sequence read " TYPE GENERIC
# OUTPUT OPTIONAL unaligned_tagged.bam
# OUTPUT OPTIONAL tagging_and_trimming_histograms.pdf 
# OUTPUT OPTIONAL preprocessed.fq.gz
# OUTPUT OPTIONAL tagging_and_trimming_summary.txt
# PARAMETER base_range_cell: "Base range for cell barcode" TYPE STRING DEFAULT 1-12 (Which bases correspond to the cell barcode)
# PARAMETER base_range_mol: "Base range for molecule barcode" TYPE STRING DEFAULT 13-20 (Which bases correspond to the molecule barcode)
# PARAMETER OPTIONAL base_quality: "Barcode quality filtering\: Minimum base quality required from the adapter bases" TYPE INTEGER DEFAULT 10 (The tool marks in the XQ tag how many bases in the barcode fall below this threshold. Barcodes and reads with bad quality bases are filtered out.)
# PARAMETER OPTIONAL sequence: "Adapter trimming\: sequence" TYPE STRING DEFAULT AAGCAGTGGTATCAACGCAGAGTGAATGGG (Adapter sequence to trim off. As a default, SMART adapter sequence.)
# PARAMETER OPTIONAL mismatches: "Adapter trimming\: Mismatches allowed in the adapter sequence" TYPE INTEGER DEFAULT 0 (How many mismatches allowed in the adapter sequence)
# PARAMETER OPTIONAL num_bases: "Adapter trimming\: Number of bases to check in adapter" TYPE INTEGER DEFAULT 5 (How many bases to check of the adapter sequence)
# PARAMETER OPTIONAL mismatches_polyA: "PolyA trimming\: Mismatches allowed in the polyA sequence" TYPE INTEGER DEFAULT 0 (How many mismatches allowed in the polyA sequence)
# PARAMETER OPTIONAL num_bases_polyA: "PolyA trimming\: Number of bases to check in polyA" TYPE INTEGER DEFAULT 6 (How many bases at least have to be A in the polyA tail)
# PARAMETER OPTIONAL minlen: "Minimum length filtering\: Minimum length of sequence reads to keep" TYPE INTEGER DEFAULT 50 (Drop the sequence read if it is below a specified length -the barcodes reads are not considered here. Trimmomatic tool is used for this step.)

# PARAMETER OPTIONAL discard_read: "Discard read" TYPE [True, False] DEFAULT False (Discard the read)
# PARAMETER OPTIONAL tag_name_cell: "Tag name for cell" TYPE [XC:XC, XM:XM] DEFAULT XC (Which BAM tag to use for the barcodes)
# OUTPUT OPTIONAL log.txt

# 2016-10-31 ML
# 2017-05-04 ML combined tools, made plots of the histograms, added Trimmomatic step
# 2017-05-18 AO Naming of outputs according to the input names NOTICE: output override

# Handle output names
# Source read_input_definitions and strip_name functions
source(file.path(chipster.common.path, "tool-utils.R"))
# read input names and strip file extension
inputnames <- read_input_definitions()
input1name <- inputnames$input.fastq.gz
input1namestripped <-strip_name(input1name)
#write the input file name into log
write(input1namestripped, file = "log.txt")

# Make a matrix of output names
# These override the default ones
outputnames <- matrix(NA, nrow=4, ncol=2)
outputnames[1,] <- c("preprocessed.fq.gz", paste(input1namestripped, ".fq.gz", sep = ""))
outputnames[2,] <- c("tagging_and_trimming_summary.txt", "tagging_and_trimming_summary.txt")
outputnames[3,] <- c("tagging_and_trimming_histograms.pdf", "tagging_and_trimming_histograms.pdf")
outputnames[4,] <- c("unaligned_tagged.bam", paste(input1namestripped, ".bam", sep =  ""))


# Write output definitions file
write_output_definitions(outputnames)



picard.binary <- file.path(chipster.tools.path, "picard-tools", "picard.jar")
path.dropseq <- c(file.path(chipster.tools.path, "drop-seq_tools"))
trimmomatic.binary <- c(file.path(chipster.tools.path, "trimmomatic", "trimmomatic-0.33.jar" ))

# STEP 1: FASTQ to BAM
# run
# single end:
if (file.exists("input2.fastq.gz")==FALSE) {
	command <- paste("java -Xmx2g -jar", picard.binary, "FastqToSam F1=input.fastq.gz O=fastq_to_bam.bam SM=something  2>> log.txt")
}

# paired end:
if (file.exists("input2.fastq.gz")==TRUE) {
	command <- paste("java -Xmx2g -jar", picard.binary, "FastqToSam F1=input.fastq.gz F2=input2.fastq.gz O=fastq_to_bam.bam SM=something  2>> log.txt")
}

#stop(paste('CHIPSTER-NOTE: ', command))
system(command)


# STEP 2: Tag BAM
# First round: cell barcode
# command start
command.start <- paste(path.dropseq, "/TagBamWithReadSequenceExtended INPUT=fastq_to_bam.bam OUTPUT=unaligned_tagged_cell.bam SUMMARY=summary_cell.txt", sep="")
# parameters
command.parameters <- paste("BASE_RANGE=", base_range_cell, " BASE_QUALITY=",base_quality ," BARCODED_READ=1 DISCARD_READ=False TAG_NAME=XC NUM_BASES_BELOW_QUALITY=1")
# run the tool
command <- paste(command.start, command.parameters, " 2>> log.txt")
system(command)
# make a plot (open the pdf)
pdf(file="tagging_and_trimming_histograms.pdf")
cell_summary <- read.table("summary_cell.txt",header = TRUE,"\t")
cell_summary2 <- data.matrix(cell_summary)
plot(cell_summary2, type = "h", col = "blue", lwd = 10, main="Number of failed cell barcodes")
# system("sed -i '1s/^/Summary_of_cell_barcodes \n/' summary_cell.txt")

# Second round: molecule barcode
# command start
command.start <- paste(path.dropseq, "/TagBamWithReadSequenceExtended INPUT=unaligned_tagged_cell.bam OUTPUT=unaligned_tagged.bam SUMMARY=summary_molecular.txt", sep="")
# parameters
command.parameters <- paste("BASE_RANGE=", base_range_mol, "BASE_QUALITY=", base_quality ," BARCODED_READ=1 DISCARD_READ=True TAG_NAME=XM NUM_BASES_BELOW_QUALITY=1")
# run the tool
command <- paste(command.start, command.parameters, " 2>> log.txt")
system(command)
# stop(paste('CHIPSTER-NOTE: ', command))
# make a plot:
molecular_summary <- read.table("summary_molecular.txt",header = TRUE,"\t")
molecular_summary2 <- data.matrix(molecular_summary)
plot(molecular_summary2, type = "h", col = "blue", lwd = 10, main="Number of failed molecule barcodes")
#system("sed -i '1s/^/Summary_of_molecular_barcodes \n/' summary_molecular.txt")

## Combine summary files:
system("sed -i \'1s/^/ Summary of cell barcodes: \\n /\' summary_cell.txt")
system("sed -i \'1s/^/ \\n \\n Summary of molecular barcodes: \\n /\' summary_molecular.txt")
system("cat summary_cell.txt summary_molecular.txt > tagging_summary.txt")


# STEP 3: Filter & trim
# FilterBAM:
command <- paste(path.dropseq, "/FilterBAM TAG_REJECT=XQ INPUT=unaligned_tagged.bam OUTPUT=unaligned_tagged_filtered.bam  2>> log.txt", sep="")
system(command)

# TrimStartingSequence:
# command start
command.start <- paste(path.dropseq, "/TrimStartingSequence INPUT=unaligned_tagged_filtered.bam OUTPUT=unaligned_tagged_trimmed.bam OUTPUT_SUMMARY=adapter_trim_summary.txt", sep="")
# parameters
command.parameters <- paste("SEQUENCE=", sequence, "MISMATCHES=", mismatches, "NUM_BASES=", num_bases)
# run the tool
command <- paste(command.start, command.parameters, " 2>> log.txt")
system(command)
# make a plot 
# pdf(file="Trimming_histogram.pdf")
polyA_trim_summary <- read.table("adapter_trim_summary.txt",header = TRUE,"\t", skip = 6)
polyA_trim_summary2 <- data.matrix(polyA_trim_summary)
plot(polyA_trim_summary2, type = "h", col = "red", lwd = 10, main="Adapter trimming")

# PolyATrimmer:
# command start
command.start <- paste(path.dropseq, "/PolyATrimmer INPUT=unaligned_tagged_trimmed.bam OUTPUT=unaligned_tagged_polyA_filtered.bam OUTPUT_SUMMARY=polyA_trimming_report.txt", sep="")
# parameters
command.parameters <- paste("MISMATCHES=", mismatches_polyA, "NUM_BASES=", num_bases_polyA)
# run the tool
command <- paste(command.start, command.parameters, " 2>> log.txt")
system(command)
# make a plot:
polyA_trim_summary <- read.table("polyA_trimming_report.txt",header = TRUE,"\t", skip = 6)
polyA_trim_summary2 <- data.matrix(polyA_trim_summary)
plot(polyA_trim_summary2, type = "h", col = "red", lwd = 5, main="polyA trimming")
dev.off() # close the pdf

# STEP 4: BAM to FASTQ
picard.binary <- file.path(chipster.tools.path, "picard-tools", "picard.jar")
# run
command <- paste("java -Xmx4g -jar", picard.binary, "SamToFastq INPUT=unaligned_tagged_polyA_filtered.bam FASTQ=preprocessed.fastq 2>> log.txt")

#stop(paste('CHIPSTER-NOTE: ', command))
system(command)


# STEP 5: Trimmomatic
# Check out if the files are compressed and if so unzip it
source(file.path(chipster.common.path, "zip-utils.R"))
unzipIfGZipFile("reads1.fastaq")
source(file.path(chipster.common.path, "zip-utils.R"))
unzipIfGZipFile("reads2.fastaq")

trim.params <- paste("")
trim.params <- paste(trim.params, "SE")
trim.params <- paste(trim.params, "-phred33")
trim.params <- paste(trim.params, "preprocessed.fastq preprocessed.fq")
step.params <- paste("")
step.params <- paste(c(step.params, " MINLEN:",  minlen), collapse="")

trimmomatic.command <- paste("java -jar", trimmomatic.binary, trim.params, step.params)
trimmomatic.command <- paste(trimmomatic.command, "1>trimlog.txt 2>> trimlog.txt")

#stop(paste('CHIPSTER-NOTE: ', trimmomatic.command))
system(trimmomatic.command)
system("gzip *.fq")

## Combine summary files:
# system("sed -i \'1s/^/ Summary of cell barcodes: \\n /\' summary_cell.txt")
system("sed -i \'1s/^/ \\n \\n Trimmomatic summary: \\n /\' trimlog.txt")
# Debug log
# system("cat log.txt > tagging_and_trimming_summary.txt")
system("cat tagging_summary.txt trimlog.txt >> tagging_and_trimming_summary.txt")

# EOF
