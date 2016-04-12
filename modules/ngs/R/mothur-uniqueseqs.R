# TOOL mothur-uniqueseqs.R: "Extract unique aligned sequences with Mothur" (Removes identical sequences from fasta-formatted alignment files. This tool is based on the Mothur package. In addition to the alignment, you need to supply the names file that was created by the tool \"Trim and filter reads with Mothur\". )
# INPUT a.fasta: "FASTA file" TYPE FASTA
# INPUT OPTIONAL a.names: "Names file" TYPE MOTHUR_NAMES
# INPUT OPTIONAL a.count_table: "Count file" TYPE GENERIC
# OUTPUT unique.fasta
# OUTPUT OPTIONAL unique.names
# OUTPUT unique.summary.tsv
# OUTPUT OPTIONAL unique.count_table


# EK 06.06.2013

# check out if the file is compressed and if so unzip it
source(file.path(chipster.common.path, "zip-utils.R"))
unzipIfGZipFile("a.fasta")

# binary
binary <- c(file.path(chipster.tools.path, "mothur", "mothur"))

# batch file
uniqueseqs.options <- ""
uniqueseqs.options <- paste(uniqueseqs.options, "unique.seqs(fasta=a.fasta")
if (file.exists("a.names")){
	uniqueseqs.options <- paste(uniqueseqs.options, " name=a.names", sep=",")
}
if (file.exists("a.count_table")){
	uniqueseqs.options <- paste(uniqueseqs.options, " count=a.count_table", sep=",")
}
uniqueseqs.options <- paste(uniqueseqs.options, ")", sep="")

# Write batch file
write(uniqueseqs.options, "batch.mth", append=F)
# command
command <- paste(binary, "batch.mth", "> log.txt 2>&1")

# run
system(command)

# Post process output
system("mv a.unique.fasta unique.fasta")
if (file.exists("a.names")){
	system("mv a.names unique.names")
}
if (file.exists("a.count_table")){
	system("mv a.count_table unique.count_table")
}

# batch file 2
write("summary.seqs(fasta=unique.fasta)", "summary.mth", append=F)

# command 2
command2 <- paste(binary, "summary.mth", "> log_raw.txt")

# run
system(command2)

# Post process output
system("grep -A 10 Start log_raw.txt > unique.summary2.tsv")
# Remove one tab to get the column naming look nice:
system("sed 's/^		/	/' unique.summary2.tsv > unique.summary.tsv")