# TOOL mothur-precluster.R: "Precluster aligned sequences" (Clusters together very similar sequences in order to remove possible sequencing errors. This tool is based on the Mothur tool pre.cluster.)
# INPUT a.fasta: "FASTA file" TYPE FASTA
# INPUT OPTIONAL a.names: "Names file" TYPE MOTHUR_NAMES
# INPUT OPTIONAL a.count_table: "Count table" TYPE GENERIC
# OUTPUT OPTIONAL preclustered.fasta
# OUTPUT OPTIONAL preclustered.names
# OUTPUT OPTIONAL preclustered-summary.tsv
# OUTPUT OPTIONAL preclustered.count_table 
# PARAMETER OPTIONAL diffs: "Number of differences allowed" TYPE INTEGER FROM 0 TO 20 DEFAULT 1 (Number of differences allowed. 1 for every 100 bp sequence is recommended)


# EK 18.06.2013
# OUTPUT OPTIONAL log.txt

# check out if the file is compressed and if so unzip it
source(file.path(chipster.common.path, "zip-utils.R"))
unzipIfGZipFile("a.fasta")

# binary
binary <- c(file.path(chipster.tools.path, "mothur", "mothur"))
# batch file
preclust.options <- ""
preclust.options <- paste(preclust.options, "pre.cluster(fasta=a.fasta")
if (file.exists("a.names")){
	preclust.options <- paste(preclust.options, " name=a.names", sep=",")
}
if (file.exists("a.count_table")){
	preclust.options <- paste(preclust.options, " count=a.count_table", sep=",")
}
preclust.options <- paste(preclust.options, ", diffs=", diffs, ")", sep="")

# Write batch file
write(preclust.options, "batch.mth", append=F)


# command
# command <- paste(binary, "batch.mth", "> log_raw.txt 2>&1")
command <- paste(binary, "batch.mth", "> log.txt")

#run
system(command)

# Post process output
system("mv a.precluster.fasta preclustered.fasta")
if (file.exists("a.names")){
	system("mv a.names preclustered.names")
}
if (file.exists("a.precluster.count_table")){
	system("mv a.precluster.count_table preclustered.count_table")
}

#stool.trim.unique.good.filter.unique.precluster.map

# batch file 2
write("summary.seqs(fasta=preclustered.fasta, count=preclustered.count_table)", "summary.mth", append=F)

# command 2
command2 <- paste(binary, "summary.mth", "> log_raw.txt")

# run
system(command2)

# Post process output
system("grep -A 10 Start log_raw.txt > preclustered-summary2.tsv")
# Remove one tab to get the column naming look nice:
system("sed 's/^		/	/' preclustered-summary2.tsv > preclustered-summary.tsv")