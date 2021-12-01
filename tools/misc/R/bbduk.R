# TOOL bbduk.R: "Bbduk example" (bbduk example)
# INPUT input.tsv: "TSV file" TYPE GENERIC
# OUTPUT output.tsv
# IMAGE comp-20.04-r-deps
# RUNTIME R-4.1.1

source(file.path(chipster.common.path, "tool-utils.R"))

runExternal("/opt/chipster/tools/bbmap/stats.sh")

runExternal(paste("mv", "input.tsv", "output.tsv"))
