# TOOL prinseq-duplicate-filter.R: "Filter reads for duplicates" (Removes duplicate reads. This tool is based on the PRINSEQ package.)
# INPUT fastqfile: "Input sequence set" TYPE GENERIC
# OUTPUT OPTIONAL accepted.fastq.gz
# OUTPUT OPTIONAL accepted.fasta.gz
# OUTPUT OPTIONAL rejected.fastq.gz
# OUTPUT OPTIONAL rejected.fasta.gz
# OUTPUT OPTIONAL filter.log
# PARAMETER derep: "Type of duplicates to filter" TYPE [1: "exact duplicate", 2: "5-prime duplicate", 3: "3-prime duplicate", 4: "reverse complement exact duplicate", 5:"reverse complement 5-prime/3-prime duplicate"] DEFAULT 1 (Type of duplicates to filter.)
# PARAMETER derep.min: "Number of allowed duplicates" TYPE INTEGER DEFAULT 2 (This option specifies the number of allowed duplicates. For example, to remove sequences that occur more than 5 times, you would specify value 6. Note that this option is used only for filtering exact duplicates or reverse complement exact duplicates.)
# PARAMETER OPTIONAL output.mode: "Results to write out" TYPE [ filt: "accepted reads only", both: "accepted and rejected reads into separate files"] DEFAULT filt (With this section you can define if the reads that get filtered out are collected to a separate file.)
# PARAMETER OPTIONAL input.mode: "Input file format" TYPE [ fq: "FASTQ", fa: "FASTA"] DEFAULT fq (Define the file format of the reads file.)
# PARAMETER OPTIONAL log.file: "Write a log file" TYPE [ n: "no", y: "yes"] DEFAULT y (Write a log file.)

# KM 17.1.2012
# AMS 17.2.2014, gzip outputs

# check out if the file is compressed and if so unzip it
source(file.path(chipster.common.lib.path, "zip-utils.R"))
unzipIfGZipFile("fastqfile")

# binary
binary.prinseq <- c(file.path(chipster.tools.path, "prinseq", "prinseq-lite.pl"))

filter.params <- paste("-derep", derep)

if (derep == "1") {
    filter.params <- paste(filter.params, "-derep_min", derep.min)
}

if (derep == "4") {
    filter.params <- paste(filter.params, "-derep_min", derep.min)
}

if (output.mode == "both") {
    filter.params <- paste(filter.params, "-out_bad rejected")
}

if (input.mode == "fq") {
    filter.command <- paste(binary.prinseq, filter.params, "-fastq fastqfile -out_good accepted -no_qual_header")
}

if (input.mode == "fa") {
    filter.command <- paste(binary.prinseq, filter.params, "-fasta fastqfile -out_good accepted")
}


if (log.file == "y") {
    system("echo Running PRINSEQ filtering with command: > filter.log")
    echo.command <- paste("echo '", filter.command, "'>> filter.log")
    system(echo.command)
    filter.command <- paste(filter.command, "-verbose 2>> filter.log")
}

system(filter.command)

system("gzip *.fastq")
system("gzip *.fasta")
