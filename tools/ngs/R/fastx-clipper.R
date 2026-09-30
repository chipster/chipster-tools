# TOOL fastx-clipper.R: "Filter reads for adapters, length and Ns with FastX" (Filters reads according to user-defined adapter sequence. Clips away the adapter,
# and optionally filters out reads that are too short after clipping or contain unknown nucleotides. Adapter-only sequences are removed in the process. This tool is based on the FASTA/Q Clipper tool of the FASTX package.)
# INPUT reads.fastq TYPE GENERIC
# OUTPUT clipped.fastq.gz
# OUTPUT clipped.log
# PARAMETER adapter: "Adapter to be removed" TYPE STRING DEFAULT CCTTAAGG (Adapter sequence that is used for filtering and that is subsequently removed.)
# PARAMETER output.options: "Output options" TYPE [clipped: "Keep only clipped reads", unclipped: "Keep only non-clipped reads", both: "Keep both clipped and non-clipped reads"] DEFAULT clipped (You can choose to keep only clipped reads (reads that contained adapter\), only non-clipped reads (reads that did not contain an adapter\) or both clipped and non-clipped reads.)
# PARAMETER minimum.alignment: "Minimum adapter alignment length" TYPE INTEGER FROM 0 DEFAULT 0 (Required minimum adapter alignment length. Maximum is adapter length. 0 means option is ignored.)
# PARAMETER OPTIONAL short: "Discard sequences shorter than" TYPE INTEGER FROM 1 DEFAULT 15 (Minimum length of sequences to keep.)
# PARAMETER OPTIONAL discard.n: "Discard sequences with Ns" TYPE [yes, no] DEFAULT yes (Keep sequences with unknown nucleotides. Default is to discard such sequences.)
# RUNTIME R-4.5.1-fastx
# TOOLS_BIN ""


# EK 27.6.2011
# AMS 01.11.2012: Added minimum.alignement and output.options options
# AMS 11.3.2014, gzip fastq outputs

# check out if the file is compressed and if so unzip it
source(file.path(chipster.common.lib.path, "zip-utils.R"))
unzipIfGZipFile("reads.fastq")

# binary
binary <- c(file.path("/opt/chipster/tools", "fastx", "bin", "fastx_clipper"))

# options
options <- paste("")
options <- paste(options, "-a", adapter)
if (minimum.alignment > 0) {
    options <- paste(options, "-M", minimum.alignment)
}
options <- paste(options, "-l", short)

if (discard.n == "no") {
    options <- paste(options, "-n")
}
if (output.options == "clipped") {
    options <- paste(options, "-c")
}
if (output.options == "unclipped") {
    options <- paste(options, "-C")
}

# command
command <- paste(binary, "-v", options, "-Q 33 -i reads.fastq -o clipped.fastq  > clipped.log")

# run
system(command)
system("gzip *.fastq")
