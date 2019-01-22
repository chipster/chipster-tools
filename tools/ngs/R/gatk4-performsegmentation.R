# TOOL gatk4-performsegmentation.R: "Segment genomic data into regions of constant copy-ratio" (GATK4)
# INPUT counts.tsv: "Tangent-normalized read counts file" TYPE GENERIC
# OUTPUT out.seg.tsv
# PARAMETER gatk4.log: "Input log2 transformed" TYPE [yes, no] DEFAULT no (Input data has had a log2 transform applied.)

source(file.path(chipster.common.path, "tool-utils.R"))
source(file.path(chipster.common.path, "zip-utils.R"))
unzipIfGZipFile("targets.tsv")

# binaries
gatk.binary <- c(file.path(chipster.tools.path, "GATK4", "gatk-launch"))

command <- paste(gatk.binary, "PerformSegmentation", "-TN counts.tsv", "-O out.seg.tsv")
if (gatk4.log == "yes"){
	command <- paste(command, "-LOG")
}

# Rscript needs to be in $PATH
set.path <-paste(sep="", "PATH=", R.home("bin"), ":$PATH")
command <- paste("bash -c '", set.path, command, "'")

runExternal(command)

# Handle output names

# read input names
inputnames <- read_input_definitions()

# Make a matrix of output names
outputnames <- matrix(NA, nrow=1, ncol=2)
outputnames[1,] <- c("out.seg.tsv", paste(strip_name(inputnames$counts.tsv), ".seg.tsv", sep=""))

# Write output definitions file
write_output_definitions(outputnames)
