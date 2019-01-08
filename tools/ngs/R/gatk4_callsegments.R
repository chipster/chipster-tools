# TOOL gatk4_callsegments.R: "Call segmented copy number variants" (Call segments as amplified, deleted, or copy number neutral given files containing tangent-normalized read counts by target and a list of segments. GATK4)
# INPUT counts.tsv: "Tangent-normalized read counts file" TYPE GENERIC
# INPUT input.seg: "Segment file" TYPE GENERIC
# OUTPUT out.called.tsv

source(file.path(chipster.common.path, "tool-utils.R"))
source(file.path(chipster.common.path, "zip-utils.R"))
unzipIfGZipFile("targets.tsv")

# binaries
gatk.binary <- c(file.path(chipster.tools.path, "GATK4", "gatk-launch"))

command <- paste(gatk.binary, "CallSegments", "-TN counts.tsv", "-S input.seg", "-O out.called.tsv")

# Rscript needs to be in $PATH
#set.path <-paste(sep="", "PATH=", R.home("bin"), ":$PATH")
#command <- paste("bash -c '", set.path, command, "'")

runExternal(command)

# Handle output names

# read input names
inputnames <- read_input_definitions()

# Make a matrix of output names
outputnames <- matrix(NA, nrow=1, ncol=2)
outputnames[1,] <- c("out.called.tsv", paste(strip_name(inputnames$counts.tsv), ".called.tsv", sep=""))

# Write output definitions file
write_output_definitions(outputnames)
