# TOOL fastx-collapser.R: "Remove duplicate reads from FASTQ" (Collapses identical sequences in a FASTQ file into a single sequence. The sequences are renamed with sequence number and the multiplicity value. This tool is based on the FASTA/Q Collapser tool of the FASTX package.)
# INPUT reads.fastq TYPE GENERIC
# OUTPUT duplicates-removed.fastq
# RUNTIME R-4.5.1-fastx
# TOOLS_BIN ""

# EK 12.1.2012

# check out if the file is compressed and if so unzip it
source(file.path(chipster.common.lib.path, "zip-utils.R"))
unzipIfGZipFile("reads.fastq")

# binary
binary <- c(file.path("/opt/chipster/tools", "fastx", "bin", "fastx_collapser"))

# command
command <- paste(binary, "-Q 33 -i reads.fastq -o duplicates-removed.fastq")

# run
system(command)
