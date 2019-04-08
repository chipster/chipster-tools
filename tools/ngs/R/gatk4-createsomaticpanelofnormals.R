# TOOL gatk4-createsomaticpanelofnormals.R: "GATK4 -Create Somatic Panel of Normals" (Create a panel of normals (PoN\) containing germline and artifactual sites for use with Mutect2. The tool takes multiple normal sample callsets produced by Mutect2's tumor-only mode and collates sites present in two or more samples into a sites-only VCF. The PoN captures common artifactual and germline variant sites. Mutect2 then uses the PoN to filter variants at the site-level. Tool is based on GATK4 createSomaticPanelOfNormals.)
# INPUT sites{...}.vcf: "VCF files" TYPE GENERIC
# OUTPUT pon.vcf.gz
# PARAMETER OPTIONAL gatk.minsamplecount: "Minumum sample count" TYPE INTEGER DEFAULT 2 (Number of samples containing a variant site required to include it in the panel of normals.) 

source(file.path(chipster.common.path, "gatk-utils.R"))
source(file.path(chipster.common.path, "tool-utils.R"))
source(file.path(chipster.common.path, "zip-utils.R"))

# binaries
gatk.binary <- c(file.path(chipster.tools.path, "GATK4", "gatk"))

command <- paste(gatk.binary, "CreateSomaticPanelOfNormals")

# Add VCF files
inputs <- paste("")
input.names <- read.table("chipster-inputs.tsv", header=F, sep="\t")
for (i in 1:nrow(input.names)) {
	if (grepl(".vcf", input.names[i,1])){
		# Index VCF
		formatGatkVcf(input.names[i,1])
		# Add input to command
		command <- paste(command, "--vcfs", paste(input.names[i,1], ".gz", sep=""))
	}	
}
command <- paste(command, "--min-sample-count", gatk.minsamplecount)
command <- paste(command, "-O pon.vcf.gz")

runExternal(command)
