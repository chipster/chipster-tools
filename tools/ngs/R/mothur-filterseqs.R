# TOOL mothur-filterseqs.R: "Filter sequence alignment" (Filters out gap columns and overhangs from a fasta formatted sequence alignment. Removing empty columns speeds up the distance calculation later. As removing columns can create new identical sequences, identical sequences are detected and removed after filtering. In addition to the FASTA file, you need to provide a count file. This tool is based on the Mothur tools filter.seqs and unique.seqs.)
# INPUT a.align: "Aligned reads in FASTA format" TYPE FASTA
# INPUT a.count_table: "Count table" TYPE MOTHUR_COUNT
# OUTPUT filtered-unique.fasta.gz
# OUTPUT filtered-log.txt
# OUTPUT filtered-unique-summary.tsv
# OUTPUT filtered-unique.count_table

# EK 05.06.2013
# ML 17.03.2017 Add optional count-table for summary file
# EK 22.03.2017 Added unique.seqs after filtering
# EK 11.05.2020 Zip output fasta

source(file.path(chipster.common.path,"tool-utils.R"))
source(file.path(chipster.common.path,"zip-utils.R"))

# check out if the file is compressed and if so unzip it
unzipIfGZipFile("a.align")

# binary
binary <- c(file.path(chipster.tools.path,"mothur","mothur"))
version <- system(paste(binary,"--version"),intern = TRUE)
documentVersion("Mothur",version)

# batch file 1
filterseqs.options <- paste("filter.seqs(fasta=a.align, vertical=T, trump=.")
filterseqs.options <- paste(filterseqs.options,", processors=",chipster.threads.max,")",sep = "")
documentCommand(filterseqs.options)
write(filterseqs.options,"batch.mth",append = FALSE)

# command
command <- paste(binary,"batch.mth","> log_raw.txt")

# run
system(command)

# result post-processing
# system("mv a.filter.fasta filtered-aligned.fasta")
system("grep -A 4 filtered log_raw.txt > filtered-log.txt")


# batch file 2
uniqueseq.options <- paste("unique.seqs(fasta=a.filter.fasta, count=a.count_table)")
documentCommand(uniqueseq.options)
write(uniqueseq.options,"batch.mth",append = FALSE)

# command 2
command2 <- paste(binary,"batch.mth","> log.txt 2>&1")

# run
system(command2)

## Post process output
system("mv a.filter.unique.fasta filtered-unique.fasta")
system("mv a.filter.count_table filtered-unique.count_table")

# batch file 3
summaryseq.options <- paste("summary.seqs(fasta=filtered-unique.fasta, count=filtered-unique.count_table")
summaryseq.options <- paste(summaryseq.options,", processors=",chipster.threads.max,")",sep = "")
documentCommand(summaryseq.options)
write(summaryseq.options,"summary.mth",append = FALSE)

# command 3
command3 <- paste(binary,"summary.mth","> log_raw.txt")

# run
system(command3)

# zip output fasta
system("gzip filtered-unique.fasta")

# Post process output
system("grep -A 10 Start log_raw.txt > filtered-summary2.tsv")
# Remove one tab to get the column naming look nice:
system("sed 's/^		/	/' filtered-summary2.tsv > filtered-unique-summary.tsv")