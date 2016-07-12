# Returns the "GATK not installed message
#
noGatkMessage <- function(){
	message <- paste("This Chipster server does not have GATK installed.\n\nGATK needs to be licensed and installed for each site running Chipster. Contact your local Chipster administrator to get GATK installed. The Chipster technical documentation has instructions on how to license, download and install GATK for Chipster.")
	return(message)
}