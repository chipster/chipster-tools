# TOOL mothur-makecontigs.R: "Combine paired reads to contigs" (Combines paired end reads to sequence contigs and puts all the resulting sequences to one fasta file. Input file is a single Tar package containing all the FASTQ files, which can be zipped. This tool is based on the Mothur tool make.contigs. Note that you can make a Tar package of your FASTQ files using the Utilities tool Make a tar package.)
# INPUT reads.tar: "Tar package containing the Fastq files" TYPE GENERIC
# OUTPUT OPTIONAL contigs.summary.tsv
# OUTPUT OPTIONAL contigs.fasta.gz
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
# AMS 16.03.2017: Changed to use single tar file as input
# OUTPUT OPTIONAL contigs.summary.tsv

source(file.path(chipster.common.path, "tool-utils.R"))
source(file.path(chipster.common.path, "zip-utils.R"))

# check out if the file is compressed and if so unzip it
unzipIfGZipFile("reads.tar")

# binary
binary <- c(file.path(chipster.tools.path, "mothur", "mothur"))
data.path <- c(file.path(chipster.tools.path, "mothur-data"))
template.path <- c(file.path(data.path, "silva.bacteria.fasta"))

# Read the contents of the tar file into a list
system("tar tf reads.tar > tar.contents")
file.list <- scan("tar.contents", what="", sep="\n")
# Check that the input is a valid tar file
if (length(file.list) == 0){
	stop(paste('CHIPSTER-NOTE: ', "It seems your input file was not a valid Tar package. Please check your input file."))
}
# Check if tar packa contains folders
if (grepl("/", file.list[1])){
	stop(paste('CHIPSTER-NOTE: ', "It seems your Tar package contains folders. The FASTQ files need to be in the root of the package, not in subfolders."))
}

# Open tar package
system("tar xf reads.tar")

# Go through the file list and gzip. Files that are already gzipped will be skipped
for (i in 1:length(file.list)) {
		system(paste("gzip", file.list[i]))
}

# Use Mothur make.file tool to generate the input list
write("make.file(inputdir=., type=gz, prefix=fastq)", "makefile.mth", append=F)
command <- paste(binary, "makefile.mth")
system(command)
system("cat *.logfile > log.tmp")


# Run Mothur make.contigs
write("make.contigs(file=fastq.files)", "makecontigs.mth", append=F)
command <- paste(binary, "makecontigs.mth", "> log2.txt")
system(command)
system("cat *.logfile >> log.tmp")

# rename the result files
system("mv fastq.trim.contigs.fasta contigs.fasta")
system("mv fastq.contigs.groups contigs.groups")

# Post process output
system("sed -n  '/Group count: / ,/Output File/p' log2.txt > log.txt")

# The summary file:
write("summary.seqs(fasta=contigs.fasta)", "summary.mth", append=F)
command <- paste(binary, "summary.mth", "> log_raw.txt")
system(command)
system("cat *.logfile >> log.tmp")

# Post process output
system("grep -A 10 Start log_raw.txt > fastq-summary2.tsv")
# Remove one tab to get the column naming look nice:
system("sed 's/^		/	/' fastq-summary2.tsv > contigs.summary.tsv.tmp")

# If contigs.summary.tsv is empty, return stderr.log
if (fileOk("contigs.summary.tsv.tmp", minlines = 1)){
	system("mv contigs.summary.tsv.tmp contigs.summary.tsv")
}else{
	system("mv log.tmp log.txt")
}

# Gzip output fasta
system("gzip contigs.fasta")