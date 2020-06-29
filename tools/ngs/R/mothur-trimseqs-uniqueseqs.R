# TOOL mothur-trimseqs-uniqueseqs.R: "Trim primers and barcodes and filter reads" (Removes primers and barcodes and filters reads for several criteria, including uniqueness. This tool is based on the Mothur tools trim.seqs, unique.seqs and count.seqs.)
# INPUT reads.fasta: "FASTA file" TYPE FASTA
# INPUT reads.oligos: "Oligos" TYPE MOTHUR_OLIGOS
# INPUT OPTIONAL reads.qual: "QUAL file" TYPE GENERIC
# OUTPUT OPTIONAL trim.unique.fasta
# OUTPUT OPTIONAL trim.groups
# OUTPUT OPTIONAL trim.unique.count_table
# OUTPUT OPTIONAL summary.trim.unique.tsv
# PARAMETER OPTIONAL flip: "Use reverse complement" TYPE [yes, no] DEFAULT no (Use reverse complement of the sequences.)
# PARAMETER OPTIONAL qaverage: "Minimum average quality of sequence" TYPE INTEGER FROM 0 TO 40 (Minimum average quality of the sequence. Sequences that have a lower average quality are dropped.)
# PARAMETER OPTIONAL qwindowaverage: "Minimum average quality of window" TYPE INTEGER (Minimum average quality score allowed over a window)
# PARAMETER OPTIONAL qwindowsize: "Window size" TYPE INTEGER (Number of bases in a window)
# PARAMETER OPTIONAL qstepsize: "Window step size" TYPE INTEGER (Number of bases to move the window over.)
# PARAMETER OPTIONAL maxambig: "Maximum ambiguous bases" TYPE INTEGER FROM 0 TO 10 (Maximum number of ambiguous bases allowed in any sequence)
# PARAMETER OPTIONAL maxhomop: "Maximum homopolymer length" TYPE INTEGER FROM 0 TO 50 (Maximum length of a homopolymer allowed in any sequence)
# PARAMETER OPTIONAL minlength: "Minimum sequence length" TYPE INTEGER FROM 0 TO 1000 (Minimum length of an allowed sequence)
# PARAMETER OPTIONAL maxlength: "Maximum sequence length" TYPE INTEGER FROM 0 TO 1000 (Maximum length of an allowed sequence)
# PARAMETER OPTIONAL pdiffs: "Maximum differences to primer sequences" TYPE INTEGER FROM 0 TO 10 (Maximum number of allowed differences to primer sequences)
# PARAMETER OPTIONAL bdiffs: "Maximum differences to barcode sequences" TYPE INTEGER FROM 0 TO 10 (Maximum number of allowed differences to barcode sequences)

# AMS 05.06.2013
# EK 11.04.2017 Renamed tool, changed names file to count file, removed qual file, modified summary to show total seqs
# OUTPUT OPTIONAL reads.trim.names
# OUTPUT OPTIONAL trim.unique.qual

source(file.path(chipster.common.path,"tool-utils.R"))
source(file.path(chipster.common.path,"zip-utils.R"))

# check out if the file is compressed and if so unzip it
unzipIfGZipFile("reads.fasta")

# binary
binary <- c(file.path(chipster.tools.path,"mothur","mothur"))
version <- system(paste(binary,"--version"),intern = TRUE)
documentVersion("Mothur",version)

# Add options
trimseqs.options <- ""
trimseqs.options <- paste(trimseqs.options,"trim.seqs(fasta=reads.fasta, oligos=reads.oligos")
if (file.exists("reads.qual")) {
  trimseqs.options <- paste(trimseqs.options," qfile=reads.qual",sep = ",")
}
if (flip == "yes") {
  trimseqs.options <- paste(trimseqs.options," flip=T",sep = ",")
}
if (!is.na(qaverage)) {
  trimseqs.options <- paste(trimseqs.options,", qaverage=",qaverage,sep = "")
}
if (!is.na(qwindowaverage)) {
  trimseqs.options <- paste(trimseqs.options,", qwindowaverage=",qwindowaverage,sep = "")
}
if (!is.na(qwindowsize)) {
  trimseqs.options <- paste(trimseqs.options,", qwindowsize=",qwindowsize,sep = "")
}
if (!is.na(qstepsize)) {
  trimseqs.options <- paste(trimseqs.options,", qstepsize=",qstepsize,sep = "")
}
if (!is.na(maxambig)) {
  trimseqs.options <- paste(trimseqs.options,", maxambig=",maxambig,sep = "")
}
if (!is.na(maxhomop)) {
  trimseqs.options <- paste(trimseqs.options,", maxhomop=",maxhomop,sep = "")
}
if (!is.na(minlength)) {
  trimseqs.options <- paste(trimseqs.options,", minlength=",minlength,sep = "")
}
if (!is.na(maxlength)) {
  trimseqs.options <- paste(trimseqs.options,", maxlength=",maxlength,sep = "")
}
if (!is.na(pdiffs)) {
  trimseqs.options <- paste(trimseqs.options,", pdiffs=",pdiffs,sep = "")
}
if (!is.na(bdiffs)) {
  trimseqs.options <- paste(trimseqs.options,", bdiffs=",bdiffs,sep = "")
}
trimseqs.options <- paste(trimseqs.options,", processors=",chipster.threads.max,sep = "")
trimseqs.options <- paste(trimseqs.options,")",sep = "")

#stop(paste('CHIPSTER-NOTE: ', trimseqs.options))

# Batch file 1 for trim.seqs and unique.seqs
documentCommand(trimseqs.options)
write(trimseqs.options,"trim.mth",append = FALSE)
uniqueseqs.options <- paste("unique.seqs(fasta=reads.trim.fasta)")
documentCommand(uniqueseqs.options)
write(uniqueseqs.options,"trim.mth",append = TRUE)
# command
command <- paste(binary,"trim.mth")
# run
system(command)
# rename fasta and groups file
system("mv reads.trim.unique.fasta trim.unique.fasta")
system("mv reads.groups trim.groups")

# Batch file 2 for count.seqs
countseqs.options <- paste("count.seqs(name=reads.trim.names, group=trim.groups)",sep = "")
documentCommand(countseqs.options)
write(countseqs.options,"batch.mth",append = FALSE)
# command
command <- paste(binary,"batch.mth",">> log.txt")
# run
system(command)
# rename count file
system("mv reads.trim.count_table trim.unique.count_table")

# Batch file 3 for summary.seqs
summaryseqs.options <- paste("summary.seqs(fasta=trim.unique.fasta, count=trim.unique.count_table")
summaryseqs.options <- paste(summaryseqs.options,", processors=",chipster.threads.max,")",sep = "")
documentCommand(summaryseqs.options)
write(summaryseqs.options,"summary.mth",append = FALSE)
# command
command <- paste(binary,"summary.mth","> log_raw.txt")
# run
system(command)

# Make trim.unique.qual
#if (file.exists("reads.trim.qual")){
#	system("grep '>' trim.unique.fasta | cut -c 2- > reads.trim.unique.list")
#	system("perl -ne 'if(/^>(\\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' reads.trim.unique.list reads.trim.qual > trim.unique.qual")
#}

# Postprocess output files
system("grep -A 10 Start log_raw.txt > summary.trim.unique2.tsv")
# Remove one tab to get the column naming look nice:
system("sed 's/^		/	/' summary.trim.unique2.tsv > summary.trim.unique.tsv")

#system("mv reads.trim.names reads.trim.names.txt")
#system("mv reads.groups reads.groups.txt")
