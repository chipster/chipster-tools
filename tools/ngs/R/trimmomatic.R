# TOOL trimmomatic.R: "Trim reads with Trimmomatic" (This tool performs a variety of trimming tasks for Illumina paired end and single end data. This tool is based on the Trimmomatic package.)
# INPUT reads1.fastaq: "Read file 1" TYPE GENERIC
# INPUT OPTIONAL reads2.fastaq: "Read file 2" TYPE GENERIC
# INPUT OPTIONAL adapters.fa: "Adapter file" TYPE GENERIC
# OUTPUT OPTIONAL trimmed.fq.gz
# OUTPUT OPTIONAL trimmed_reads1_paired.fq.gz
# OUTPUT OPTIONAL trimmed_reads1_unpaired.fq.gz
# OUTPUT OPTIONAL trimmed_reads2_paired.fq.gz
# OUTPUT OPTIONAL trimmed_reads2_unpaired.fq.gz
# OUTPUT OPTIONAL trimlog.txt
# PARAMETER OPTIONAL adapter.file: "Adapter set" TYPE [none: "none", TruSeq2-SE.fa: "TruSeq2-SE", TruSeq3-SE.fa: "TruSeq3-SE", TruSeq2-PE.fa: "TruSeq2-PE", TruSeq3-PE.fa: "TruSeq3-PE", TruSeq3-PE-2.fa: "TruSeq3-PE-2", NexteraPE-PE.fa: "NexteraPE-PE"] DEFAULT none (Cut adapter and other Illumina-specific sequences from the read. You will also need to provide cutting parameters.)
# PARAMETER OPTIONAL illuminaclip: "Adapter clipping parameters" TYPE STRING (You will need to supply minimum three parameters: <seed mismatches>:<palindrome clip threshold>:<simple clip threshold>. There are also two optional parameters that affect palindrome mode only: <min adapter length> and <keep both reads>. Value for <keep both reads> is given as true/false. Values are separated by colons, e.g. 2:30:10 or 2:30:10:1:true)
# PARAMETER OPTIONAL phred.scale: "Quality scale used in the fastq file" TYPE [phred33: "phred + 33", phred64: "phred + 64"] DEFAULT phred33 (Quality scale used in the fastq file.)
# PARAMETER OPTIONAL leading: "Minimum quality to keep a leading base" TYPE INTEGER (Cut bases off the start of a read, if below a threshold quality.)
# PARAMETER OPTIONAL trailing: "Minumum quality to keep a trailing base" TYPE INTEGER (Cut bases off the end of a read, if below a threshold quality.)
# PARAMETER OPTIONAL crop: "Number of bases to keep from the start" TYPE INTEGER (Cut the read to a specified length.)
# PARAMETER OPTIONAL headcrop: "Number of bases to remove from the start" TYPE INTEGER (Cut the specified number of bases from the start of the read.)
# PARAMETER OPTIONAL slidingwindow: "Sliding window trimming parameters" TYPE STRING (Perform a sliding window trimming from the 5' end, cutting once the average quality within the window falls below a threshold. Required parameter are <window size>:<required quality>. Values are separated by colons e.g. 4:15)
# PARAMETER OPTIONAL maxinfo: "Adaptive quality trimming parameters" TYPE STRING (An adaptive quality trimmer which balances the benefits of retaining longer reads against the costs of retaining bases with errors. Two parameters need to be provided: <target length>:<strictness>. Strictness is a decimal value between 0 and 1. Values are separated by colons, e.g. 36:0.8)
# PARAMETER OPTIONAL avgqual: "Minimum average quality of reads to keep" TYPE INTEGER (Drop the read if the average quality is below the specified level.)
# PARAMETER OPTIONAL minlen: "Minimum length of reads to keep" TYPE INTEGER (Drop the read if it is below a specified length.)
# PARAMETER OPTIONAL unpaired: "Write unpaired reads in a separate file" TYPE [yes, no] DEFAULT no (When using Trimmomatic for paired reads, some reads might lose their mate. These unpaired reads are not included in the result file, but you can request to get them in a separate file if needed.)
# PARAMETER OPTIONAL logfile: "Write a log file" TYPE [yes, no] DEFAULT no (Write a log file.)
# STORAGE 1000

# AMS 2014.04.08
# MK, 2014.12.05, corrected typo: avqual => avgqual. Corrected bug in initialisation of adapter.file parameter
# AMS 2014.11.27, corrected bug: trimmomatic was always run in SE mode
# ML, 2015.12.17, added option to use own adapter files
# EK, 2024.01.15, added storage

source(file.path(chipster.common.lib.path, "tool-utils.R"))

# Check out if the files are compressed and if so unzip it
source(file.path(chipster.common.lib.path, "zip-utils.R"))
unzipIfGZipFile("reads1.fastaq")
source(file.path(chipster.common.lib.path, "zip-utils.R"))
unzipIfGZipFile("reads2.fastaq")


trimmomatic.binary <- c(file.path(chipster.tools.path, "trimmomatic", "trimmomatic-0.33.jar"))

# Parameters
trim.params <- paste("")
if (file.exists("reads2.fastaq")) {
    # Paired end reads
    trim.params <- paste(trim.params, "PE")
    trim.params <- paste(trim.params, "-threads", chipster.threads.max)
    if (phred.scale == "phred33") {
        trim.params <- paste(trim.params, "-phred33")
    } else {
        trim.params <- paste(trim.params, "-phred64")
    }
    trim.params <- paste(trim.params, "reads1.fastaq reads2.fastaq trimmed_reads1_paired.fq trimmed_reads1_unpaired.fq trimmed_reads2_paired.fq trimmed_reads2_unpaired.fq")
} else {
    # Single end reads
    trim.params <- paste(trim.params, "SE")
    trim.params <- paste(trim.params, "-threads", chipster.threads.max)
    if (phred.scale == "phred33") {
        trim.params <- paste(trim.params, "-phred33")
    } else {
        trim.params <- paste(trim.params, "-phred64")
    }
    trim.params <- paste(trim.params, "reads1.fastaq trimmed.fq")
}

# Trimming steps

step.params <- paste("")
if (adapter.file != "none") {
    adapter.path <- c(file.path(chipster.tools.path, "trimmomatic", "adapters", adapter.file))
    if (!nchar(illuminaclip) > 0) {
        stop("CHIPSTER-NOTE: You need to provide the required parameters for the adapter clipping to work.")
    }
    step.params <- paste(c(step.params, " ILLUMINACLIP:", adapter.path, ":", illuminaclip), collapse = "")
}
if (file.exists("adapters.fa")) {
    if (!nchar(illuminaclip) > 0) {
        stop("CHIPSTER-NOTE: You need to provide the required parameters for the adapter clipping to work.")
    }
    if (adapter.file != "none") {
        stop("CHIPSTER-NOTE: Choose either one of the adapter sets or use your own adapter file, don't do both.")
    }
    step.params <- paste(c(step.params, " ILLUMINACLIP:", "adapters.fa", ":", illuminaclip), collapse = "")
}
if (!is.na(leading)) {
    step.params <- paste(c(step.params, " LEADING:", leading), collapse = "")
}
if (!is.na(trailing)) {
    step.params <- paste(c(step.params, " TRAILING:", trailing), collapse = "")
}
if (!is.na(crop)) {
    step.params <- paste(c(step.params, " CROP:", crop), collapse = "")
}
if (!is.na(headcrop)) {
    step.params <- paste(c(step.params, " HEADCROP:", headcrop), collapse = "")
}
if (nchar(slidingwindow) > 0) {
    step.params <- paste(c(step.params, " SLIDINGWINDOW:", slidingwindow), collapse = "")
}
if (nchar(maxinfo) > 0) {
    step.params <- paste(c(step.params, " MAXINFO:", maxinfo), collapse = "")
}
if (!is.na(avgqual)) {
    step.params <- paste(c(step.params, " AVGQUAL:", avgqual), collapse = "")
}
if (!is.na(minlen)) {
    step.params <- paste(c(step.params, " MINLEN:", minlen), collapse = "")
}

# Check that at leaat one
if (!nchar(step.params) > 0) {
    stop("CHIPSTER-NOTE: No trimming or clipping steps selected. Select at least one.")
}
trimmomatic.command <- paste("java -jar", trimmomatic.binary, trim.params, step.params)
if (logfile == "yes") {
    trimmomatic.command <- paste(trimmomatic.command, "1>trimlog.txt 2>> trimlog.txt")
}
# stop(paste('CHIPSTER-NOTE: ', trimmomatic.command))

documentCommand(trimmomatic.command)
system(trimmomatic.command)

system("gzip *.fq")


# Handle output names
source(file.path(chipster.common.lib.path, "tool-utils.R"))

# read input names
inputnames <- read_input_definitions()

# Make a matrix of output names
outputnames <- matrix(NA, nrow = 6, ncol = 2)

base1 <- strip_name(inputnames$reads1.fastaq)
if (file.exists("reads2.fastaq")) {
    base2 <- strip_name(inputnames$reads2.fastaq)
} else {
    base2 <- ""
}

outputnames[1, ] <- c("trimmed.fq.gz", paste(base1, "_trimmed.fq.gz", sep = ""))
outputnames[2, ] <- c("trimmed_reads1_paired.fq.gz", paste(base1, "_trimmed.fq.gz", sep = ""))
outputnames[3, ] <- c("trimmed_reads1_unpaired.fq.gz", paste(base1, "_unpaired_trimmed.fq.gz", sep = ""))
outputnames[4, ] <- c("trimmed_reads2_paired.fq.gz", paste(base2, "_trimmed.fq.gz", sep = ""))
outputnames[5, ] <- c("trimmed_reads2_unpaired.fq.gz", paste(base2, "_unpaired_trimmed.fq.gz", sep = ""))
outputnames[6, ] <- c("trimlog.txt", paste(base1, "_trimlog.txt", sep = ""))

# Delete files containing unpaired reads if they are not needed
if (unpaired == "no") {
    system("rm *unpaired.fq.gz")
}

# Write output definitions file
write_output_definitions(outputnames)
