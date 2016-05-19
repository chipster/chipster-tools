# TOOL mothur-countseqs.R: "Count sequences with Mothur" (Counts the number of sequences represented by the representative sequence in a name file. If a group file is given, it will also provide the group count breakdown. This tool is based on the Mothur package.)
# INPUT a.names: "Names file" TYPE MOTHUR_NAMES
# INPUT a.fasta: "FASTA file" TYPE FASTA
# INPUT OPTIONAL a.groups: "Groups file" TYPE MOTHUR_GROUPS
# OUTPUT counts.count_table
# OUTPUT count-summary.tsv

# ML 30.03.2016

# binary
binary <- c(file.path(chipster.tools.path, "mothur", "mothur"))

# batch file
if (file.exists("a.groups")){
	write(paste("count.seqs(name=a.names, group=a.groups)", sep=""), "batch.mth", append=F)
}else {
	write(paste("count.seqs(name=a.names)", sep=""), "batch.mth", append=F)
}

# command
command <- paste(binary, "batch.mth", "> log.txt")

#run
system(command)

# Post process output
system("mv a.count_table counts.count_table")


# batch file 2
write("summary.seqs(fasta=a.fasta, count=counts.count_table)", "summary.mth", append=F)

# command 2
command2 <- paste(binary, "summary.mth", "> log_raw.txt")

# run
system(command2)

# Post process output
system("grep -A 10 Start log_raw.txt > count-summary2.tsv")
# Remove one tab to get the column naming look nice:
system("sed 's/^		/	/' count-summary2.tsv > count-summary.tsv")