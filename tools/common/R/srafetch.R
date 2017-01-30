# TOOL srafetch.R: "Retrieve FASTQ files from SRA database" (Retrieve reads in FASTQ format from the SRA database based on the entry ID or name.)
# OUTPUT OPTIONAL srafetch.log
# OUTPUT sra_reads_{...}.fastqc.gz: "FASTQ files"
# PARAMETER entry_id: "Name or SRR ID of the SRA entry" TYPE STRING DEFAULT "entry" (Give the SRR id of the SRA entry to be retrieved. For example: SRR000021) 
# PARAMETER dump: "Sequences to dump" TYPE [all: "All", aligned: "Only aligned sequences", unaligned: "Only unaligned sequences"  ] DEFAULT all (Define the reads to be retrieved from the SRA entry)

# KM 08.11.2014
# ML 19.12.2016
# TH 20.1.2017

dump.param <- ""
if (dump == "aligned") {
	dump.param <- "--aligned"
} else if (dump == "unaligned") {
	dump.param <- "--unaligned"
}

sra.path <- file.path(chipster.tools.path, "sratoolkit", "bin")
sra.binary <- file.path(sra.path, "fastq-dump")
command.full <- paste(sra.binary, dump.param, '--split-files --gzip', entry_id,  '1>>srafetch.log 2>>srafetch.log')
cat(command.full, "\n", file="srafetch.log", append=TRUE)
system(command.full)

# check that we got something
fastq.files <- Sys.glob("*.fastq.gz")
if (length(fastq.files) == 0) {
	system('echo srafetch.log:')
	system('cat srafetch.log')
	stop("Fetching fastq files failed, see the end of the details for more information")
}

# rename outputs and create output names file
source(file.path(chipster.common.path, "tool-utils.R"))
outputnames <- matrix(NA, nrow=length(fastq.files), ncol=2)
for (i in 1:length(fastq.files)) {
	original.name <- fastq.files[i]
	new.name <- paste("sra_reads_", i, ".fastqc.gz", sep ="")
	outputnames[i,] <- c(new.name, original.name)
	system(paste("mv", original.name, new.name))
}
write_output_definitions(outputnames)

