# TOOL fastqc_multiqc_storage.R: "Read quality with MultiQC for input data bigger than 200 GB" (The tool runs FastQC on multiple FASTQ files, and then combines the reports using MultiQC. Input can be gzipped FASTQ files or tar files containing FASTQ files. The total size of input files can be bigger than 200 GB. Please make sure you don't have duplicate FASTQ file names. Run the tool once for all samples, not separately for each file. This tool is based on the FastQC and MultiQC packages.)
# INPUT reads{...}.fq: "FASTQ files" TYPE GENERIC
# OUTPUT OPTIONAL multiqc_report.html
# OUTPUT OPTIONAL error_log.txt
# RUNTIME R-4.1.1-fastqc
# STORAGE 1000

# 2018.09.10 AMS
# 2022.12.02 ES updated to new version

source(file.path(chipster.common.lib.path, "tool-utils.R"))
# source(file.path(chipster.common.lib.path, "zip-utils.R"))

# Rename input files to get the original banes in the report
input.names <- read.table("chipster-inputs.tsv", header = FALSE, sep = "\t")

system("cat chipster-inputs.tsv")
system("ls -l")

for (i in 1:nrow(input.names)) {
  system(paste("mv --backup=numbered", input.names[i, 1], input.names[i, 2]))
}
system("ls -l")

# Try opening the input files as tar packages.
system("ls *.tar *.tar.* *.tgz |xargs -i tar xf {} --xform='s#^.+/##x' 2>/dev/null")

# binary
fastqc.binary <- file.path(chipster.tools.path, "fastqc", "fastqc")
multiqc.binary <- file.path(chipster.tools.path, "multiqc", "multiqc")

# Document versions
version <- system(paste(fastqc.binary, "--version"), intern = TRUE)
documentVersion("FastQC", version)
version <- system(paste(multiqc.binary, "--version"), intern = TRUE)
documentVersion("MultiQC", version)

fastq.files <- system("ls *.fastq *.fastq.gz *.fq *.fq.qz", intern = TRUE)

# commands
fastqc.command <- paste(fastqc.binary, "-f fastq --noextract ", paste(fastq.files, collapse = " "))
multiqc.command <- paste(multiqc.binary, "--module fastqc .")

documentCommand(fastqc.command)
documentCommand(multiqc.command)

# run
runExternal(fastqc.command)
runExternal(multiqc.command)



if (fileNotOk("multiqc_report.html")) {
  system("mv stderr.log error_log.txt")
}
