# TOOL mothur-uniqueseqs.R: "Extract unique sequences" (Identifies identical sequences and keep only one representative sequence and the number of times this sequence occures in each sample. If names file is given, the tool also counts the number of sequences represented by the representative sequence in a names file, generating a counts_table. If a groups file is also given, it will also provide the group count breakdown. This tool is based on the Mothur tools unique.seqs and count.seqs.)
# INPUT a.fasta: "FASTA file" TYPE FASTA
# INPUT OPTIONAL a.names: "Names file" TYPE MOTHUR_NAMES
# INPUT OPTIONAL a.count_table: "Count file" TYPE MOTHUR_COUNT
# INPUT OPTIONAL a.groups: "Groups file" TYPE MOTHUR_GROUPS
# OUTPUT unique.fasta
# OUTPUT unique.summary.tsv
# OUTPUT OPTIONAL unique.count_table


# OUTPUT OPTIONAL unique.names

# EK 06.06.2013
# ML 17.3.2017 Combine with countseqs.R 

# check out if the file is compressed and if so unzip it
source(file.path(chipster.common.path, "zip-utils.R"))
unzipIfGZipFile("a.fasta")

# binary
binary <- c(file.path(chipster.tools.path, "mothur", "mothur"))

# batch file 1
uniqueseqs.options <- ""
uniqueseqs.options <- paste(uniqueseqs.options, "unique.seqs(fasta=a.fasta")
if (file.exists("a.names")){
	uniqueseqs.options <- paste(uniqueseqs.options, " name=a.names", sep=",")
}
if (file.exists("a.count_table")){
	uniqueseqs.options <- paste(uniqueseqs.options, " count=a.count_table", sep=",")
}
uniqueseqs.options <- paste(uniqueseqs.options, ")", sep="")

# Write batch file 1
write(uniqueseqs.options, "batch.mth", append=F)
# command
command <- paste(binary, "batch.mth", "> log.txt 2>&1")
# run
system(command)

## Post process output
system("mv a.unique.fasta unique.fasta")
#if (file.exists("a.names")){
#	system("mv a.names unique.names")
#}


# batch file 2 count.seqs, if names and/or group file exists
if (file.exists("a.names")){
	if (file.exists("a.groups")){
	write(paste("count.seqs(name=a.names, group=a.groups)", sep=""), "batch.mth", append=F)
	}else {
	write(paste("count.seqs(name=a.names)", sep=""), "batch.mth", append=F)
	}
	
# command 2
command <- paste(binary, "batch.mth", "> log.txt")
#run
system(command)
}

# Post process output
if (file.exists("a.count_table")){
	system("mv a.count_table unique.count_table")
}

# batch file 3
if (file.exists("unique.count_table")){
	write("summary.seqs(fasta=unique.fasta, count=unique.count_table)", "summary.mth", append=F)
} else {
		write("summary.seqs(fasta=unique.fasta)", "summary.mth", append=F)
		}

# command 3
command2 <- paste(binary, "summary.mth", "> log_raw.txt")
# run
system(command2)

# Post process output
system("grep -A 10 Start log_raw.txt > unique.summary2.tsv")
# Remove one tab to get the column naming look nice:
system("sed 's/^		/	/' unique.summary2.tsv > unique.summary.tsv")