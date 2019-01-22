# TOOL gatk4-plotsegmentedcopyratio.R: "Create plots of denoised and segmented copy ratio" (GATK4)
# INPUT tn.tsv: "Tangent-normalized coverage" TYPE GENERIC
# INPUT ptn.tsv: "Pre-tangent-normalized coverage" TYPE GENERIC
# INPUT input.seg: "Segment file" TYPE GENERIC
# INPUT dict: "Sequence dictionary file" TYPE GENERIC
# OUTPUT OPTIONAL chipster_Before_After_CR_Lim_4.png  
# OUTPUT OPTIONAL gatk4_dQc.txt
# OUTPUT OPTIONAL gatk4_postQc.txt
# OUTPUT OPTIONAL gatk4_scaled_dQc.txt
# OUTPUT OPTIONAL gatk4_Before_After.png
# OUTPUT OPTIONAL gatk4_FullGenome.png
# OUTPUT OPTIONAL gatk4_preQc.tx
# PARAMETER gatk4.log: "Input log2 transformed" TYPE [yes, no] DEFAULT no (Input data has had a log2 transform applied.)

source(file.path(chipster.common.path, "tool-utils.R"))
source(file.path(chipster.common.path, "zip-utils.R"))

unzipIfGZipFile("tn.tsv")
unzipIfGZipFile("ptn.tsv")

# binaries
gatk.binary <- c(file.path(chipster.tools.path, "GATK4", "gatk-launch"))

command <- paste(gatk.binary, "PlotSegmentedCopyRatio", "-TN tn.tsv", "-PTN ptn.tsv", "-S input.seg", "-O plots", "-pre gatk4", "-SD dict")
if (gatk4.log == "yes"){
	command <- paste(command, "-LOG")
}

# Rscript needs to be in $PATH
set.path <-paste(sep="", "PATH=", R.home("bin"), ":$PATH")
command <- paste("bash -c '", set.path, command, "'")

# Output folder must exist
system ("mkdir plots")

Sys.setenv("DISPLAY"=":0")

runExternal(command)

system("cp plots/* .")

# Handle output names

# read input names
#inputnames <- read_input_definitions()

# Make a matrix of output names
#outputnames <- matrix(NA, nrow=1, ncol=2)
#outputnames[1,] <- c("out.seg.tsv", paste(strip_name(inputnames$counts.tsv), ".seg.tsv", sep=""))

# Write output definitions file
#write_output_definitions(outputnames)
