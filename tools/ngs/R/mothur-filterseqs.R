# TOOL mothur-filterseqs.R: "Filter sequence alignment" (Filters out gap columns and overhangs from a fasta formatted sequence alignment. Removing empty columns speeds up the distance calculation later. As removing columns can create new identical sequences, identical sequences are detected and removed after filtering. In addition to the FASTA file, you need to provide a count file. This tool is based on the Mothur tools filter.seqs and unique.seqs.)
# INPUT a.align: "Aligned reads in FASTA format" TYPE GENERIC
# INPUT a.count_table: "Count table" TYPE MOTHUR_COUNT
# OUTPUT filtered-unique.fasta
# OUTPUT filtered-log.txt
# OUTPUT filtered-unique-summary.tsv
# OUTPUT filtered-unique.count_table

# EK 05.06.2013
# ML 17.03.2017 Add optional count-table for summary file
# EK 22.03.2017 Added unique.seqs after filtering

# check out if the file is compressed and if so unzip it
source(file.path(chipster.common.path, "zip-utils.R"))
unzipIfGZipFile("a.align")

# binary
binary <- c(file.path(chipster.tools.path, "mothur", "mothur"))

# batch file 1
write("filter.seqs(fasta=a.align, vertical=T, trump=.)", "batch.mth", append=F)

# command
command <- paste(binary, "batch.mth", "> log_raw.txt")

# run
system(command)

# result post-processing
# system("mv a.filter.fasta filtered-aligned.fasta")
system("grep -A 4 filtered log_raw.txt > filtered-log.txt")


# batch file 2
write("unique.seqs(fasta=a.filter.fasta, count=a.count_table)", "batch.mth", append=F)

# command 2
command2 <- paste(binary, "batch.mth", "> log.txt 2>&1")

# run
system(command2)

## Post process output
system("mv a.filter.unique.fasta filtered-unique.fasta")
system("mv a.filter.count_table filtered-unique.count_table")

# batch file 3

write("summary.seqs(fasta=filtered-unique.fasta, count=filtered-unique.count_table)", "summary.mth", append=F)

# command 3
command3 <- paste(binary, "summary.mth", "> log_raw.txt")

# run
system(command3)

# Post process output
system("grep -A 10 Start log_raw.txt > filtered-summary2.tsv")
# Remove one tab to get the column naming look nice:
system("sed 's/^		/	/' filtered-summary2.tsv > filtered-unique-summary.tsv")