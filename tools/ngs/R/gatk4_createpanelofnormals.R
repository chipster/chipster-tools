# TOOL gatk4_createpanelofnormals.R: "Create Panel of Normals \(PoN\)" (GATK4)
# INPUT counts.tsv: "proportional read counts as a TSV file" TYPE GENERIC
# OUTPUT pon.pon
# OUTPUT OPTIONAL removed_samples.txt
# OUTPUT OPTIONAL target_weights.txt

source(file.path(chipster.common.path, "tool-utils.R"))
source(file.path(chipster.common.path, "zip-utils.R"))

unzipIfGZipFile("counts.tsv")

# binaries
gatk.binary <- c(file.path(chipster.tools.path, "GATK4", "gatk-launch"))

command <- paste(gatk.binary, "CreatePanelOfNormals", "-I counts.tsv", "-O pon.pon", "--disableSpark")

runExternal(command)

# Handle output names

# read input names
inputnames <- read_input_definitions()

# Make a matrix of output names
outputnames <- matrix(NA, nrow=1, ncol=2)
outputnames[1,] <- c("pon.pon", paste(strip_name(inputnames$counts.tsv), ".pon", sep=""))

# Write output definitions file
write_output_definitions(outputnames)

if (fileOk("pon.pon.removed_samples.txt", minsize=1 )){
	system("mv pon.pon.removed_samples.txt removed_samples.txt")
}
if (fileOk("pon.pon.target_weights.txt", minsize=1 )){
	system("mv pon.pon.target_weights.txt target_weights.txt")
}
