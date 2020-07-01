# TOOL mothur-summaryseqs.R: "Summarize sequences" (Given a fasta file with unaligned or aligned sequences, provides summary statistics on sequence start and end coordinates, length, number of ambiguous bases, and homopolymer length. This tool is based on the Mothur tool summary.seqs.)
# INPUT reads.fasta: "FASTA file" TYPE FASTA
# OUTPUT summary.tsv

# AMS 4.6.2013
# EK 27.6.2013 Changes to description and output

source(file.path(chipster.common.path,"tool-utils.R"))
source(file.path(chipster.common.path,"zip-utils.R"))

# check out if the file is compressed and if so unzip it
unzipIfGZipFile("reads.fasta")

# binary
binary <- c(file.path(chipster.tools.path,"mothur","mothur"))
version <- system(paste(binary,"--version"),intern = TRUE)
documentVersion("Mothur",version)

# batch file
summaryseqs.options <- paste("summary.seqs(fasta=reads.fasta, processors=",chipster.threads.max,")",sep = "")
documentCommand(summaryseqs.options)
write(summaryseqs.options,"batch.mth",append = FALSE)

# command
command <- paste(binary,"batch.mth","> log_raw.txt")

# run
system(command)

# Post process output
# system("mv reads.summary summary.tsv")
# Choose the summary lines:
system("grep -A 10 Start log_raw.txt > summary_2.tsv")
# Remove one tab to get the column naming look nice:
system("sed 's/^		/	/' summary_2.tsv > summary.tsv")
