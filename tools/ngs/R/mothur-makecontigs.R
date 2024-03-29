# TOOL mothur-makecontigs.R: "Combine paired reads to contigs with Mothur" (Combines paired reads to sequence contigs within each sample, and puts all the resulting sequences in one fasta file. Input file is a single Tar package containing all the FASTQ files, which can be gzipped. You can make a Tar package of your FASTQ files using the Utilities tool Make a tar package. The tool tries to assign the FASTQ files into samples based on the file names, but you can also provide a file containing this information, please see the manual. We recommend you to use tool "Combine paired reads to contigs with VSEARCH" because it performs better)
# INPUT reads.tar: "Tar package containing the FASTQ files" TYPE GENERIC
# INPUT OPTIONAL input_list: "List of FASTQ files by sample" TYPE GENERIC
# OUTPUT contigs.summary.tsv
# OUTPUT contigs.fasta.gz
# OUTPUT contigs.count_table
# OUTPUT contig.numbers.txt
# OUTPUT samples.fastqs.txt
# RUNTIME R-4.1.1

# ML 02.03.2016
# AMS 16.03.2017: Changed to use single tar file as input
# ES 1.12.2022 Changed to use new mothur version 1.48, produce count file not a gorup file

# OUTPUT log3.txt

source(file.path(chipster.common.lib.path, "tool-utils.R"))
source(file.path(chipster.common.lib.path, "zip-utils.R"))

# check out if the file is compressed and if so unzip it
unzipIfGZipFile("reads.tar")

# binary
binary <- c(file.path(chipster.tools.path, "mothur", "mothur"))
data.path <- c(file.path(chipster.tools.path, "mothur-data"))
template.path <- c(file.path(data.path, "silva.bacteria.fasta"))
version <- system(paste(binary, "--version"), intern = TRUE)
documentVersion("Mothur", version)

# Read the contents of the tar file into a list
system("tar tf reads.tar > tar.contents")
file.list <- scan("tar.contents", what = "", sep = "\n")
# Check that the input is a valid tar file
if (length(file.list) == 0) {
  stop(paste("CHIPSTER-NOTE: ", "It seems your input file is not a valid Tar package. Please check your input file."))
}
# Check if tar packa contains folders
if (grepl("/", file.list[1])) {
  stop(paste("CHIPSTER-NOTE: ", "It seems your Tar package contains folders. The FASTQ files need to be in the root of the package, not in subfolders."))
}

# Open tar package
system("tar xf reads.tar")

# Go through the file list and gzip. Files that are already gzipped will be skipped
for (i in 1:length(file.list)) {
  system(paste("gzip", file.list[i]))
}

# Use input_list if provided. Else use Mothur make.file tool to generate the input list
if (fileOk("input_list")) {
  system("mv input_list fastq.files")
} else {
  write("make.file(inputdir=., type=gz, prefix=fastq)", "makefile.mth", append = FALSE)
  command <- paste(binary, "makefile.mth")
  system(command)
  system("cat *.logfile > log.tmp")
  # Remove full paths from the Mothur-generated input list to make it more readable
  system("for line in $( cat fastq.files ); do echo `basename $line`; done | paste - - - > samples.fastqs.txt")
}

# Run Mothur make.contigs
makecontigs.options <- paste("make.contigs(file=fastq.files")
makecontigs.options <- paste(makecontigs.options, ", processors=", chipster.threads.max, ")", sep = "")
documentCommand(makecontigs.options)
write(makecontigs.options, "makecontigs.mth", append = FALSE)
command <- paste(binary, "makecontigs.mth", "> log2.txt")
system(command)
system("cat *.logfile >> log.tmp")

# rename the result files
system("mv fastq.trim.contigs.fasta contigs.fasta")


# Run Mothur count.seqs
countseqs.options <- paste("count.seqs(count=fastq.contigs.count_table, compress=f)", sep = "")
documentCommand(countseqs.options)
write(countseqs.options, "countseqs.mth", append = FALSE)
command <- paste(binary, "countseqs.mth", "> log3.txt")
system(command)
system("cat *.logfile >> log.tmp")
system("mv fastq.contigs.full.count_table contigs.count_table")

# Post process output
system("sed -n  '/Group count: / ,/Output File/p' log2.txt > contig.numbers.txt")

# The summary file:
summaryseqs.options <- paste("summary.seqs(fasta=contigs.fasta")
summaryseqs.options <- paste(summaryseqs.options, ", processors=", chipster.threads.max, ")", sep = "")
documentCommand(summaryseqs.options)
write(summaryseqs.options, "summary.mth", append = FALSE)
command <- paste(binary, "summary.mth", "> log_raw.txt")
system(command)
system("cat *.logfile >> log.tmp")

# Post process output
system("grep -A 10 Start log_raw.txt > fastq-summary2.tsv")
# Remove one tab to get the column naming look nice:
system("sed 's/^		/	/' fastq-summary2.tsv > contigs.summary.tsv.tmp")

# If contigs.summary.tsv is empty, return the log instead
if (fileOk("contigs.summary.tsv.tmp", minlines = 1)) {
  system("mv contigs.summary.tsv.tmp contigs.summary.tsv")
} else {
  system("mv log.tmp log.txt")
}

# Gzip output fasta
system("gzip contigs.fasta")
