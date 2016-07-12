# TOOL no_gatk_example.R: "No GATK example" (Example)
# INPUT reference.fasta: "Reference genome" TYPE GENERIC


# AMS 11.07.2016

# binary
gatk.binary <- c(file.path(chipster.tools.path, "GATK", " GenomeAnalysisTK-eioo.jar"))

# Check if GATK is installed
if (file.exists(gatk.binary) == FALSE){
	source(file.path(chipster.common.path, "gatk-utils.R"))
	message <- noGatkMessage()
	stop(paste('CHIPSTER-NOTE: ', message))
}
