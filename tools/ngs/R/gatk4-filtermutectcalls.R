# TOOL gatk4-filtermutectcalls.R: "GATK4 -Filter Mutect2 calls" (Filter variants in a Mutect2 VCF callset. Tool is based on GATK4 FilterMutectCalls.)
# INPUT input.vcf: "Mutect2 calls" TYPE GENERIC
# INPUT contaminationtable: "Contamination table" TYPE GENERIC
# OUTPUT OPTIONAL mutect2_filtered.vcf
# OUTPUT OPTIONAL statistics.tsv
# OUTPUT OPTIONAL gatk_log.txt

source(file.path(chipster.common.path, "gatk-utils.R"))
source(file.path(chipster.common.path, "tool-utils.R"))
source(file.path(chipster.common.path, "zip-utils.R"))

unzipIfGZipFile("contaminationtable")

# binaries
gatk.binary <- c(file.path(chipster.tools.path, "GATK4", "gatk"))

# Preprocess VCF
formatGatkVcf("input.vcf")

# Run FilterMutectCalls
command <- paste(gatk.binary, "FilterMutectCalls", "-V input.vcf.gz", "-contamination-table contaminationtable", "--stats statistics.tsv", "-O mutect2_filtered.vcf")
runExternal(command)

# Return error message if no result
if (fileNotOk("mutect2_filtered.vcf")){
	system("mv stderr.log gatk_log.txt")
}

# read input names
inputnames <- read_input_definitions()

basename <- strip_name(inputnames$input.vcf)

# Make a matrix of output names
outputnames <- matrix(NA, nrow=2, ncol=2)
outputnames[1,] <- c("mutect2_filtered.vcf", paste(basename, "_filtered.vcf", sep=""))
outputnames[2,] <- c("statistics.tsv", paste(basename, "_statistics.tsv", sep=""))

# Write output definitions file
write_output_definitions(outputnames)