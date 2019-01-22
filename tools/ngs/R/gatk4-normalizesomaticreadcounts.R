# TOOL gatk4-normalizesomaticreadcounts.R: "Normalize PCOV read counts using a panel of normals" (GATK4)
# INPUT counts.tsv: "Read counts as a TSV file" TYPE GENERIC
# INPUT pon.pon: "PoN file" TYPE GENERIC
# OUTPUT ptn.tsv
# OUTPUT tn.tsv

source(file.path(chipster.common.path, "tool-utils.R"))
source(file.path(chipster.common.path, "zip-utils.R"))

unzipIfGZipFile("counts.tsv")

# binaries
gatk.binary <- c(file.path(chipster.tools.path, "GATK4", "gatk-launch"))

command <- paste(gatk.binary, "NormalizeSomaticReadCounts", "-I counts.tsv", "-PON pon.pon", "-PTN ptn.t.tsv", "-TN tn.t.tsv")

runExternal(command)

# Cleand headers to make compatible with Chipster tsv viewer
system("grep -v '#' ptn.t.tsv > ptn.tsv")
system("grep -v '#' tn.t.tsv > tn.tsv")

# Handle output names

# read input names
inputnames <- read_input_definitions()

# Make a matrix of output names
outputnames <- matrix(NA, nrow=2, ncol=2)
outputnames[1,] <- c("ptn.tsv", paste(strip_name(inputnames$counts.tsv), ".ptn.tsv", sep=""))
outputnames[2,] <- c("tn.tsv", paste(strip_name(inputnames$counts.tsv), ".tn.tsv", sep=""))

# Write output definitions file
write_output_definitions(outputnames)
