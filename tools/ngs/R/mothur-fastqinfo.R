# TOOL mothur-fastqinfo.R: "Split FASTQ file to FASTA and QUAL files" (Given a FASTQ file, produces a FASTA file and QUAL file. This tool is based on the Mothur tool fastq.info.)
# INPUT a.fastq: "FASTQ file" TYPE GENERIC
# OUTPUT reads.fasta
# OUTPUT basequalities.qual
# OUTPUT reads.summary.tsv
# RUNTIME R-4.1.1

# EK 11.04.2017
# ES 28.12.2022 use new version of Mothur

source(file.path(chipster.common.lib.path, "tool-utils.R"))
source(file.path(chipster.common.lib.path, "zip-utils.R"))

# check out if the file is compressed and if so unzip it
unzipIfGZipFile("a.fastq")

# binary
# binary <- c(file.path(chipster.tools.path,"mothur-1.44.3","mothur"))
binary <- c(file.path(chipster.tools.path, "mothur", "mothur"))
version <- system(paste(binary, "--version"), intern = TRUE)
documentVersion("Mothur", version)

# batch file
fastqinfo.options <- paste("fastq.info(fastq=a.fastq)")
documentCommand(fastqinfo.options)
write(fastqinfo.options, "batch.mth", append = FALSE)

# command
command <- paste(binary, "batch.mth", "> log_raw.txt")

# run
system(command)

# rename output files
system("mv a.fasta reads.fasta")
system("mv a.qual basequalities.qual")

# batch file 2
summaryseqs.options <- paste("summary.seqs(fasta=reads.fasta, processors=", chipster.threads.max, ")", sep = "")
documentCommand(summaryseqs.options)
write(summaryseqs.options, "batch.mth", append = FALSE)

# command
command <- paste(binary, "batch.mth", "> log_raw.txt")

# run
system(command)

# Post process output
# Choose the summary lines:
system("grep -A 10 Start log_raw.txt > summary_2.tsv")
# Remove one tab to get the column naming look nice:
system("sed 's/^		/	/' summary_2.tsv > reads.summary.tsv")
