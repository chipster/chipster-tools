# TOOL gatk4_padtargets.R: "Pad target intervals" (GATK4)
# INPUT targets.tsv: "Targets TSV" TYPE GENERIC
# OUTPUT padded.tsv

source(file.path(chipster.common.path, "tool-utils.R"))
source(file.path(chipster.common.path, "zip-utils.R"))
unzipIfGZipFile("targets.tsv")

# binaries
gatk.binary <- c(file.path(chipster.tools.path, "GATK4", "gatk-launch"))

command <- paste(gatk.binary, "PadTargets", "-T targets.tsv", "-O padded.tsv")

runExternal(command)

# Handle output names

# read input names
inputnames <- read_input_definitions()

# Make a matrix of output names
outputnames <- matrix(NA, nrow=1, ncol=2)
outputnames[1,] <- c("padded.tsv", paste("pad_", inputnames$targets.tsv, sep=""))

# Write output definitions file
write_output_definitions(outputnames)
