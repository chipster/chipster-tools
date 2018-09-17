# TOOL fastqc_multiqc.R: "Read quality with MultiQC for many FASTQ files" (The tool runs FastQC on multiple FASTQ files, and then combines the reports using MultiQC. Input file is a single Tar package containing all the FASTQ files, which can be gzipped. This tool is based on the FastQC and MultiQC packages.)
# INPUT reads.tar: "Tar package containing FASTQ files" TYPE GENERIC
# OUTPUT OPTIONAL multiqc_report.html
# OUTPUT OPTIONAL error_log.txt

# 2018.09.10 AMS

source(file.path(chipster.common.path, "tool-utils.R"))


# Read the contents of the tar file into a list
system("tar tf reads.tar > tar.contents 2>> log.txt")
file.list <- scan("tar.contents", what="", sep="\n")

# Check that the input is a valid tar file
if (length(file.list) == 0){
	stop(paste('CHIPSTER-NOTE: ', "It seems your input file is not a valid Tar package. Please check your input file."))
}

# Open tar package. Folders are flattened.
system("tar xf reads.tar --xform='s#^.+/##x' 2>> log.txt")	

# binary
fastqc.binary <- file.path(chipster.tools.path, "FastQC", "fastqc")
# Needs a link in tools like RSeQC
multiqc.binary <- file.path(chipster.tools.path, "Python-2.7.12", "bin", "multiqc")

# commands
fastqc.command <- paste(fastqc.binary, "-f fastq --noextract *")
multiqc.command <- paste(multiqc.binary, ".")

# run
runExternal(fastqc.command)
runExternal(multiqc.command)

if(fileNotOk("multiqc_report.html")){
	system("mv stderr.log error_log.txt")
}

