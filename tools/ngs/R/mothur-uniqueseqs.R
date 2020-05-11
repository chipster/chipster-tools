# TOOL mothur-uniqueseqs.R: "Extract unique sequences" (Given a fasta file and group file, this tool identifies identical sequences and keeps only one representative sequence. It stores the number of times each representative sequence occurs in each sample in a count_table file. This tool is based on the Mothur tools unique.seqs and count.seqs.)
# INPUT a.fasta: "FASTA file" TYPE FASTA
# INPUT a.groups: "Groups file" TYPE MOTHUR_GROUPS
# OUTPUT unique.fasta.gz
# OUTPUT OPTIONAL unique.summary.tsv
# OUTPUT OPTIONAL unique.count_table


# OUTPUT OPTIONAL log.txt

# EK 06.06.2013
# ML 17.03.2017 Combine with countseqs.R 
# ML 22.03.2017 Remove extra input options -these functionalities are moved to the trimming tool
# EK 11.05.2020 Zip output fasta

# check out if the file is compressed and if so unzip it
source(file.path(chipster.common.path, "zip-utils.R"))
unzipIfGZipFile("a.fasta")

# binary
binary <- c(file.path(chipster.tools.path, "mothur", "mothur"))

# batch file 1
uniqueseqs.options <- ""
uniqueseqs.options <- paste(uniqueseqs.options, "unique.seqs(fasta=a.fasta")
uniqueseqs.options <- paste(uniqueseqs.options, ")", sep="")

# Write batch file 1
write(uniqueseqs.options, "batch.mth", append=F)
# command
command <- paste(binary, "batch.mth", "> log.txt 2>&1")
# run
system(command)

## Post process output
system("mv a.unique.fasta unique.fasta")



# batch file 2 count.seqs -creates a new count_table
write(paste("count.seqs(name=a.names, group=a.groups)", sep=""), "batch.mth", append=F)
	
# command 2
command <- paste(binary, "batch.mth", ">> log.txt")
#run
system(command)
# Output File Names: 
#		a.count_table
	
# Post process output -note, same name as with the input file, thus inside the "if"!
if (file.exists("a.count_table")){
	system("mv a.count_table unique.count_table")
}	


# batch file 3 -summary
write("summary.seqs(fasta=unique.fasta, count=unique.count_table)", "summary.mth", append=F)


# command 3
command2 <- paste(binary, "summary.mth", "> log_raw.txt")

# run
system(command2)

# zip output fasta
system("gzip unique.fasta")

# Post process output
system("grep -A 10 Start log_raw.txt > unique.summary2.tsv")
# Remove one tab to get the column naming look nice:
system("sed 's/^		/	/' unique.summary2.tsv > unique.summary.tsv")