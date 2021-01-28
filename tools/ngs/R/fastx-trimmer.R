# TOOL fastx-trimmer.R: "Trim reads with FastX" (Trims reads to a user-specified length. This tool is based on the FASTA/Q Trimmer tool of the FASTX package.)
# INPUT reads.fastq TYPE GENERIC
# OUTPUT trimmed.fastq.gz 
# PARAMETER first: "First base to keep" TYPE INTEGER FROM 1 TO 100 DEFAULT 1 (First base to keep.)
# PARAMETER last: "Last base to keep" TYPE INTEGER FROM 1 TO 300 DEFAULT 75 (Last base to keep.)
# PARAMETER quality.format: "Quality value format used" TYPE [sanger: Sanger, illuminaold: "Illumina GA v1.3-1.5"] DEFAULT sanger (What quality encoding is used in your FASTQ file. Select Sanger if your data comes from Illumina 1.8 or later, SOLiD or 454.)




# EK 17.6.2011
# AMS 11.3.2014, gzip fastq outputs

source(file.path(chipster.common.path,"tool-utils.R"))
source(file.path(chipster.common.path, "zip-utils.R"))

# check out if the file is compressed and if so unzip it
unzipIfGZipFile("reads.fastq")

# binary
binary <- c(file.path(chipster.tools.path, "fastx", "bin", "fastx_trimmer"))
version <- system(paste(binary,"-h | sed -n 2p"),intern = TRUE)
documentVersion("fastx_trimmer",version)

# command
quality.scale <- ifelse(quality.format == "sanger", "-Q 33", "")
command <- paste(binary, "-f", first, "-l", last, quality.scale, "-i reads.fastq -o trimmed.fastq")
documentCommand(command)

# run
system(command)
system("gzip *.fastq")

# Determine base name
inputnames <- read_input_definitions()
basename <- strip_name(inputnames$reads.fastq)

# Make a matrix of output names
outputnames <- matrix(NA,nrow = 1,ncol = 2)
outputnames[1,] <- c("trimmed.fastq.gz",paste(basename,".trimmed.fq.gz",sep = ""))

# Write output definitions file
write_output_definitions(outputnames)