# TOOL gatk4-filtermutectcalls.R: "Filter Mutect2 calls with GATK4" (Uses FilterMutectCalls.)
# INPUT input.vcf: "Mutect calls" TYPE GENERIC
# INPUT contaminationtable: "Contamination table" TYPE GENERIC
# OUTPUT OPTIONAL mutect2_filtered.vcf
# OUTPUT OPTIONAL gatk_log.txt

source(file.path(chipster.common.path, "tool-utils.R"))
source(file.path(chipster.common.path, "zip-utils.R"))

unzipIfGZipFile("reference")

# binaries
gatk.binary <- c(file.path(chipster.tools.path, "GATK4", "gatk-launch"))

# Run FilterMutectCalls
command <- paste(gatk.binary, "FilterMutectCalls", "-O mutect2_filtered.vcf", "-V input.vcf", "-contaminationTable contaminationtable")
runExternal(command)


# Return error message if no result
if (fileNotOk("mutect2_filtered.vcf")){
	system("mv stderr.log gatk_log.txt")
}

