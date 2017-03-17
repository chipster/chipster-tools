# TOOL mothur-makecontigs.R: "Combine paired reads to contigs" (Combines paired end reads to sequence contigs and puts all the resulting sequences to one fasta file. You need a list of your samples in .files -format as an input. This tool is based on the Mothur tool make.contigs.)
# INPUT stability.files: "stability.files" TYPE MOTHUR_STABILITY
# INPUT sample{...}.fastq: "Fastq files" TYPE GENERIC
# OUTPUT contigs.summary.tsv
# OUTPUT OPTIONAL contigs.fasta
# OUTPUT OPTIONAL contigs.groups
# OUTPUT OPTIONAL log.txt


#Output File Names: 
#inputs.trim.contigs.fasta
#inputs.trim.contigs.qual
#inputs.contigs.report
#inputs.scrap.contigs.fasta
#inputs.scrap.contigs.qual
#inputs.contigs.groups
# OUTPUT log_raw.txt
# OUTPUT OPTIONAL fastq.contigs.report
# OUTPUT OPTIONAL fastq.trim.contigs.qual


# ML 02.03.2016
# OUTPUT OPTIONAL contigs.summary.tsv



# Use the chipster-inputs.tsv to match the chipster input names to actual file names:
input.names <- read.table("chipster-inputs.tsv", header=F, sep="\t")
orig.file.names <- read.table("stability.files", header=F, sep="\t")

# Sanity check the input data:
if (nrow(orig.file.names)*2+1 > nrow(input.names)) {
	stop("CHIPSTER-NOTE: You have more files listed in your .files table than you have as inputs!")
}
if (nrow(orig.file.names)*2+1 < nrow(input.names)) {
	stop("CHIPSTER-NOTE: You have some extra files as inputs!")
}
# Create a new stability.files using the input file names:
new.files <- mat.or.vec(nrow(orig.file.names), 3)
for (i in 1:nrow(orig.file.names)) {
	new.files[i,1] <- as.character(orig.file.names[i,1])
	new.files[i,2] <- as.character(input.names[!is.na(pmatch(input.names[,2], orig.file.names[i,2])),1])
	new.files[i,3] <- as.character(input.names[!is.na(pmatch(input.names[,2], orig.file.names[i,3])),1])
}


# write.table(input.names, "inputs.tsv", col.names=T, row.names=T, sep="\t", quote=F)
write.table(new.files, "fastq.files", col.names=F, row.names=F, sep="\t", quote=F)


# binary
binary <- c(file.path(chipster.tools.path, "mothur", "mothur"))
data.path <- c(file.path(chipster.tools.path, "mothur-data"))
template.path <- c(file.path(data.path, "silva.bacteria.fasta"))

# batch file
write("make.contigs(file=fastq.files)", "batch.mth", append=F)

# command
command <- paste(binary, "batch.mth", "> log2.txt")

# run
system(command)

# rename the result files
system("mv fastq.trim.contigs.fasta contigs.fasta")
system("mv fastq.contigs.groups contigs.groups")


# Post process output
system("sed -n  '/Group count: / ,/Output File/p' log2.txt > log.txt")
#$ sed -n '/WORD1/,/WORD2/p' /path/to/file
#$ sed -n '/FOO/,/BAR/p' test.txt

# The summary file:
write("summary.seqs(fasta=contigs.fasta)", "summary.mth", append=F)

# command
command2 <- paste(binary, "summary.mth", "> log_raw.txt")

# run
system(command2)

# Post process output
system("grep -A 10 Start log_raw.txt > fastq-summary2.tsv")
# Remove one tab to get the column naming look nice:
system("sed 's/^		/	/' fastq-summary2.tsv > contigs.summary.tsv")
