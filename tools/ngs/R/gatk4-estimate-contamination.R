# TOOL gatk4-estimate-contamination.R: "GATK4 -Estimate contamination" (Calculates the fraction of reads coming from cross-sample contamination, given results from GetPileupSummaries. Input is a pileup table in TSV format. You can also provide an optional matched control pileup table. The resulting contamination table is used with FilterMutectCalls. Tool is based on GATK4 CalculateContamination.)
# INPUT tumor.tsv: "Tumor pileup table" TYPE GENERIC
# INPUT OPTIONAL normal.tsv: "Normal pileup table" TYPE GENERIC
# OUTPUT OPTIONAL CalculateContamination.tsv
# OUTPUT OPTIONAL gatk_log.txt

source(file.path(chipster.common.path, "tool-utils.R"))
source(file.path(chipster.common.path, "zip-utils.R"))

unzipIfGZipFile("tumor.tsv")
unzipIfGZipFile("normal.tsv")

# binaries
gatk.binary <- c(file.path(chipster.tools.path, "GATK4", "gatk"))

options <-""

options <- paste("-I tumor.tsv")

if (fileOk("normal.tsv")){
	options <- paste(options, "-matched normal.tsv")	
}

# Run CalculateContamination
command <- paste(gatk.binary, "CalculateContamination", options, "-O CalculateContamination.tsv")
runExternal(command)

# Return error message if no result
if (fileNotOk("CalculateContamination.tsv")){
	system("mv stderr.log gatk_log.txt")
}

# read input names
inputnames <- read_input_definitions()

basename <- strip_name(inputnames$tumor.tsv)
basename <- remove_postfix(basename, "_getpileupsummaries")

# Make a matrix of output names
outputnames <- matrix(NA, nrow=1, ncol=2)
if (fileOk("normal.tsv")){
	outputnames[1,] <- c("CalculateContamination.tsv", paste(basename, "_pair_calculatecontamination.tsv", sep=""))
}else{
	outputnames[1,] <- c("CalculateContamination.tsv", paste(basename, "_calculatecontamination.tsv", sep=""))
}
# Write output definitions file
write_output_definitions(outputnames)