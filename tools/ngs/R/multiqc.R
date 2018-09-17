# TOOL multiqc.R: "Combine reports using MultiQC" (This tool combines reports from supported tools such as FastQC. This tool is based on the MultiQC package.)
# INPUT reports{...}: "Reports" TYPE GENERIC
# OUTPUT OPTIONAL multiqc_report.html
# OUTPUT OPTIONAL error_log.txt

# 2018.09.10 AMS

source(file.path(chipster.common.path, "tool-utils.R"))

# All inputs are expected to be .tgz
input.names <- read.table("chipster-inputs.tsv", header=F, sep="\t")
for (i in 1:nrow(input.names)) {
	system(paste("tar xf",(input.names[i,1])))	
}


# binary
# Needs a link in tools like RSeQC
binary <- file.path(chipster.tools.path, "Python-2.7.12", "bin", "multiqc")

# command
command <- paste(binary, ".")

# run
runExternal(command)

if(fileNotOk("multiqc_report.html")){
	system("mv stderr.log error_log.txt")
}

