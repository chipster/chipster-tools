# TOOL fastx-quality-filter.R: "Filter reads for quality with FastX" (Filters out reads which contain quality scores below the user-specified criteria. This tool is based on the FASTQ Quality Filter tool of the FASTX package.)
# INPUT reads.fastq TYPE GENERIC
# OUTPUT quality-filtered.fastq.gz
# OUTPUT quality-filtered.log
# PARAMETER quality: "Quality cut-off value" TYPE INTEGER FROM 1 TO 100 DEFAULT 20 (What is the minimum quality score to keep.)
# PARAMETER percentage: "Minimum percent of bases that must have that quality" TYPE INTEGER FROM 1 TO 100 DEFAULT 90 (Percent of bases in sequence that must have quality equal to or higher than the cut-off value.)
# RUNTIME R-4.5.1-fastx
# TOOLS_BIN ""

# EK 28.6.2011
# AMS 11.3.2014, gzip fastq outputs

# check out if the file is compressed and if so unzip it
source(file.path(chipster.common.lib.path, "zip-utils.R"))
unzipIfGZipFile("reads.fastq")

# binary
binary <- c(file.path("/opt/chipster/tools", "fastx", "bin", "fastq_quality_filter"))

# command
command <- paste(binary, "-v", "-q", quality, "-p", percentage, "-Q 33 -i reads.fastq -o quality-filtered.fastq > quality-filtered.log")

# run
system(command)
system("gzip *.fastq")
