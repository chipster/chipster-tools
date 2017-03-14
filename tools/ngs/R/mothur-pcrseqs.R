# TOOL mothur-pcrseqs.R: "Extract amplified region from reference alignment" (Extract your amplified region from the Silva reference alignment. This tool is based on the Mothur tool pcr.seqs.)
# INPUT silva.bacteria.fasta: "reference FASTA file" TYPE FASTA
# OUTPUT custom.reference.fasta
# OUTPUT custom.reference.summary.tsv
# PARAMETER OPTIONAL start: "Start" TYPE INTEGER (Start point of your region of interest)
# PARAMETER OPTIONAL end: "End" TYPE INTEGER (End point of your region of interest)
# PARAMETER OPTIONAL keepdots: "Remove leading and trailing dots" TYPE [yes, no] DEFAULT yes (Remove leading and trailing dots.)

# ML 23.3.2016

## check out if the file is compressed and if so unzip it
#source(file.path(chipster.common.path, "zip-utils.R"))
#unzipIfGZipFile("reads.fasta")

# binary
binary <- c(file.path(chipster.tools.path, "mothur", "mothur"))
data.path <- c(file.path(chipster.tools.path, "mothur-data"))
template.path <- c(file.path(data.path, "silva.bacteria.fasta"))

# Add options
pcrseqs.options <- ""
pcrseqs.options <- paste(pcrseqs.options, "pcr.seqs(fasta=silva.bacteria.fasta", sep="")
if (!is.na(start)){
	pcrseqs.options <- paste(pcrseqs.options, ", start=", start, sep="")
}
if (!is.na(end)){
	pcrseqs.options <- paste(pcrseqs.options, ", end=", end, sep="")
}
if (keepdots=="yes"){
	pcrseqs.options <- paste(pcrseqs.options, ", keepdots=F", sep="")
}
if (keepdots=="no"){
	pcrseqs.options <- paste(pcrseqs.options, ", keepdots=T", sep="")
}
pcrseqs.options <- paste(pcrseqs.options, ")", sep="")

# Write batch file
write(pcrseqs.options, "pcrseq.mth", append=F)

# command
command <- paste(binary, "pcrseq.mth")

# run
system(command)

# rename the file
system("mv silva.bacteria.pcr.fasta custom.reference.fasta")

# batch file 2
write("summary.seqs(fasta=custom.reference.fasta)", "summary.mth", append=F)

# command
command2 <- paste(binary, "summary.mth", "> log_raw.txt")

# run
system(command2)

# Post process output
system("grep -A 10 Start log_raw.txt > custom.reference.summary2.tsv")
# Remove one tab to get the column naming look nice:
system("sed 's/^		/	/' custom.reference.summary2.tsv > custom.reference.summary.tsv")