# TOOL bbduk.R: "Bbduk example" (bbduk example)
# INPUT input.tsv: "TSV file" TYPE GENERIC
# OUTPUT output.tsv
# RUNTIME R-4.1.1

source(file.path(chipster.common.path, "tool-utils.R"))

bbmap_dir <- paste(chipster.tools.path, "bbmap", sep="/")
stats_path <- paste(bbmap_dir, "stats.sh", sep="/")

# run with bash...
runExternal(paste("bash", "-c", stats_path))

# ... because the script doesn't work without it:
# runExternal(stats_path)

runExternal(paste("mv", "input.tsv", "output.tsv"))
