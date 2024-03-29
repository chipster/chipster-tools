# TOOL mothur-precluster.R: "Precluster aligned sequences" (Clusters together very similar sequences in order to remove possible sequencing errors. In addition to the fasta file, you need to give count_table file or names file. Note that sample names should not contain hyphens. This tool is based on the Mothur tool pre.cluster.)
# INPUT a.fasta: "FASTA file" TYPE FASTA
# INPUT OPTIONAL a.names: "Names file" TYPE MOTHUR_NAMES
# INPUT OPTIONAL a.count_table: "Count table" TYPE MOTHUR_COUNT
# OUTPUT preclustered.fasta.gz
# OUTPUT OPTIONAL preclustered.names
# OUTPUT preclustered-summary.tsv
# OUTPUT OPTIONAL preclustered.count_table
# PARAMETER OPTIONAL diffs: "Number of differences allowed" TYPE INTEGER FROM 0 TO 20 DEFAULT 1 (Number of differences allowed for every 100 bases. 1 for every 100 bp sequence is recommended.)
# RUNTIME R-4.1.1
# SLOTS 4

# EK 18.06.2013
# EK 06.05.2020 Type for input count_table changed
# EK 11.05.2020 Zip output fasta
# ES 1.12.2022 Changed to use new mothur version 1.48

# OUTPUT OPTIONAL log.txt
# add if for names file in summary

source(file.path(chipster.common.lib.path, "tool-utils.R"))
source(file.path(chipster.common.lib.path, "zip-utils.R"))

# check out if the file is compressed and if so unzip it
unzipIfGZipFile("a.fasta")

# binary
binary <- c(file.path(chipster.tools.path, "mothur", "mothur"))
# binary <- c(file.path(chipster.tools.path,"mothur-1.44.3","mothur"))
version <- system(paste(binary, "--version"), intern = TRUE)
documentVersion("Mothur", version)

# batch file
preclust.options <- ""
preclust.options <- paste(preclust.options, "pre.cluster(fasta=a.fasta")
if (file.exists("a.names")) {
  preclust.options <- paste(preclust.options, " name=a.names", sep = ",")
}
if (file.exists("a.count_table")) {
  preclust.options <- paste(preclust.options, " count=a.count_table", sep = ",")
}
preclust.options <- paste(preclust.options, ", diffs=", diffs, sep = "")
preclust.options <- paste(preclust.options, ", processors=", chipster.threads.max, ")", sep = "")

# Write batch file
documentCommand(preclust.options)
write(preclust.options, "batch.mth", append = FALSE)


# command
# command <- paste(binary, "batch.mth", "> log_raw.txt 2>&1")
command <- paste(binary, "batch.mth", "> log.txt")

# run
system(command)

# Post process output
system("mv a.precluster.fasta preclustered.fasta")
if (file.exists("a.names")) {
  system("mv a.names preclustered.names")
}
if (file.exists("a.precluster.count_table")) {
  system("mv a.precluster.count_table preclustered.count_table")
}

# stool.trim.unique.good.filter.unique.precluster.map

# batch file 2
# write("summary.seqs(fasta=preclustered.fasta, count=preclustered.count_table)", "summary.mth", append=F)
summaryseqs.options <- paste("summary.seqs(fasta=preclustered.fasta")
if (file.exists("preclustered.count_table")) {
  summaryseqs.options <- paste(summaryseqs.options, ", count=preclustered.count_table", sep = "")
}
summaryseqs.options <- paste(summaryseqs.options, ", processors=", chipster.threads.max, ")", sep = "")
documentCommand(summaryseqs.options)
write("summary.seqs(fasta=preclustered.fasta)", "summary.mth", append = FALSE)

# command 2
command2 <- paste(binary, "summary.mth", "> log_raw.txt")

# run
system(command2)

# zip output fasta
system("gzip preclustered.fasta")

# Post process output
system("grep -A 10 Start log_raw.txt > preclustered-summary2.tsv")
# Remove one tab to get the column naming look nice:
system("sed 's/^		/	/' preclustered-summary2.tsv > preclustered-summary.tsv")
