# TOOL bwa-mem-with-index-building.R: "BWA MEM for single or paired end reads and own genome" (This tool uses BWA MEM to align single-end or paired-end reads to a provided reference genome. You need to supply the single end or paired-end reads in FASTQ format. If you provide one reads file as input, single-end alignment is performed. If you provide paired-end reads files as input, paired-end alignment is performed. You also need to provide the reference genome as input either in FASTA format or as a tar package with a BWA MEM index.)
# INPUT reads{...}.fq: "Reads" TYPE FASTQ
# INPUT OPTIONAL genome.txt: "Reference genome" TYPE GENERIC
# INPUT OPTIONAL reads1.txt: "List of read 1 files" TYPE TEXT_STRICT
# INPUT OPTIONAL reads2.txt: "List of read 2 files" TYPE TEXT_STRICT
# OUTPUT bwa.bam
# OUTPUT bwa.log
# OUTPUT OPTIONAL bwa.bam.bai
# OUTPUT OPTIONAL bwa_index.tar
# PARAMETER OPTIONAL index.file: "Create index file" TYPE [index_file: "Create index file", no_index: "No index file"] DEFAULT no_index (Creates index file for BAM. By default no index file.)
# PARAMETER OPTIONAL minseedlen: "Minimum seed length" TYPE INTEGER DEFAULT 19 (Matches shorter than this will be missed when looking for maximal exact matches or MEMs in the first alignment phase.)
# PARAMETER OPTIONAL bandwith: "Maximum gap length" TYPE INTEGER DEFAULT 100 (Gaps longer than this will not be found. Note also scoring matrix and hit length affect the maximum gap length, in addition to this band width parameter.)
# PARAMETER OPTIONAL matchscore: "Match score" TYPE INTEGER DEFAULT 1 (Score for a matching base.)
# PARAMETER OPTIONAL mismatchscore: "Mismatch penalty" TYPE INTEGER DEFAULT 4 (Penalty for a mismatching base.)
# PARAMETER OPTIONAL gapopen: "Gap opening penalty" TYPE INTEGER DEFAULT 6 (Gap opening penalty.)
# PARAMETER OPTIONAL gapextension: "Gap extension penalty" TYPE INTEGER DEFAULT 1 (Gap extension penalty.)
# PARAMETER OPTIONAL clippenalty: "Penalty for end clipping" TYPE INTEGER DEFAULT 5 (Penalty for 5\'- and 3\'-end clipping. When performing the Smith-Waterman extension of the seed alignments, BWA-MEM keeps track of the best score reaching the end of the read. If this score is larger than the best SW score minus the clipping penalty, clipping will not be applied.)
# PARAMETER OPTIONAL rgid: "Read group identifier" TYPE STRING (Read group identifier. If you want to add the read group line in the BAM file, you have to give this information.)
# PARAMETER OPTIONAL rgsm: "Sample name for read group" TYPE STRING (The name of the sample sequenced in this read group. Note that you have to fill in also the read group identifier parameter for the read group information to appear in the BAM file.)
# PARAMETER OPTIONAL rgpl: "Platform for read group" TYPE [ none: "Not defined", ILLUMINA, SOLID, LS454, HELICOS, PACBIO] DEFAULT none (Platform\/technology used to produce the read. Note that you have to fill in also the read group identifier parameter for the read group information to appear in the BAM file.)
# PARAMETER OPTIONAL rglb: "Library identifier for read group" TYPE STRING (DNA preparation library identifier. The Mark Duplicates tool uses this field to determine which read groups might contain molecular duplicates, in case the same DNA library was sequenced on multiple lanes. Note that you have to fill in also the read group identifier parameter for the read group information to appear in the BAM file.)
# RUNTIME R-4.1.1

# KM 1.9.2015

source(file.path(chipster.common.lib.path, "tool-utils.R"))
source(file.path(chipster.common.lib.path, "bam-utils.R"))
source(file.path(chipster.common.lib.path, "zip-utils.R"))

# check out if the file is compressed and if so unzip it
input.names <- read.table("chipster-inputs.tsv", header = F, sep = "\t")
for (i in 1:nrow(input.names)) {
    unzipIfGZipFile(input.names[i, 1])
}



# read input names
inputnames <- read_input_definitions()

# bwa binary
bwa.binary <- file.path(chipster.tools.path, "bwa", "bwa")
# bwa mem binary
bwa.mem.binary <- paste(bwa.binary, "mem")
bwa.index.binary <- file.path(chipster.module.path, "shell", "check_bwa_index.sh")
samtools.binary <- c(file.path(chipster.tools.path, "samtools", "bin", "samtools"))

genome.filetype <- system("file -b genome.txt | cut -d ' ' -f2", intern = TRUE)
hg_ifn <- ("")
echo.command <- paste("echo Host genome file type", genome.filetype, " > bwa.log")
runExternal(echo.command)

new_index_created <- ("no")
# case 1. Ready calculated indexes in tar format
if (genome.filetype == "tar") {
    runExternal("echo Extracting tar formatted gemome index file >> bwa.log")
    runExternal("tar -tf genome.txt >> bwa.log")
    check.command <- paste(bwa.index.binary, "genome.txt | tail -1 ")
    bwa.genome <- system(check.command, intern = TRUE)
    if (bwa.genome == "wrong_tar_content") {
        stop("CHIPSTER-NOTE: The selected genome file does not contain BWA indexes.")
    }
    runExternal("ls -l >> bwa.log")
    # case 2. Fasta file
} else {
    # check sequece file type
    emboss.path <- file.path(chipster.tools.path, "emboss-20.04", "bin")
    options(scipen = 999)
    inputfile.to.check <- ("genome.txt")
    sfcheck.binary <- file.path(chipster.module.path, "../misc/shell/sfcheck.sh")
    sfcheck.command <- paste(sfcheck.binary, emboss.path, inputfile.to.check)
    str.filetype <- system(sfcheck.command, intern = TRUE)

    if (str.filetype == "Not an EMBOSS compatible sequence file") {
        stop("CHIPSTER-NOTE: Your reference genome is not a sequence file that is compatible with the tool you try to use")
    }
    # Do indexing
    runExternal("echo Calculating gemome indexes >> bwa.log")
    check.command <- paste(bwa.index.binary, "genome.txt -tar| tail -1 ")
    bwa.genome <- system(check.command, intern = TRUE)
    cp.command <- paste("cp ", bwa.genome, "_bwa_index.tar ./bwa_index.tar ", sep = "")
    runExternal(cp.command)
    new_index_created <- ("yes")
}
echo.command <- paste("echo Internal genome name:", bwa.genome, " >> bwa.log")
runExternal(echo.command)

# bwa
command.start <- (bwa.mem.binary)

bwa.parameters <- paste("-M", "-k", minseedlen, "-w", bandwith, "-A", matchscore, "-B", mismatchscore, "-O", gapopen, "-E", gapextension, "-L", clippenalty)

# Read group definition
if (nchar(rgid) > 0) {
    rg.string <- ("'@RG")
    rg.string <- paste(rg.string, "\\tID:", rgid, sep = "")
    if (nchar(rgsm) > 0) {
        rg.string <- paste(rg.string, "\\tSM:", rgsm, sep = "")
    }
    if (rgpl != "none") {
        rg.string <- paste(rg.string, "\\tPL:", rgpl, sep = "")
    }
    if (nchar(rglb) > 0) {
        rg.string <- paste(rg.string, "\\tLB:", rglb, sep = "")
    }
    rg.string <- paste(rg.string, "'", sep = "")
    bwa.parameters <- paste(bwa.parameters, "-R", rg.string)
}

# Read group error message
if (nchar(rgid) < 1) {
    if (nchar(rgsm) > 0) {
        stop("CHIPSTER-NOTE: Please define identifier for read group")
    }
    if (rgpl != "none") {
        stop("CHIPSTER-NOTE: Please define identifier for read group")
    }
    if (nchar(rglb) > 0) {
        stop("CHIPSTER-NOTE: Please define identifier for read group")
    }
}

# Input files
if (fileOk("reads1.txt")) {
    if (fileNotOk("reads2.txt")) {
        # Case: One list file -> multiple single-end alignments
        reads1.list <- make_input_list("reads1.txt")
        paired.end <- FALSE
    } else {
        # Case: Two list files -> multiple paired-end alignments
        reads1.list <- make_input_list("reads1.txt")
        reads2.list <- make_input_list("reads2.txt")
        paired.end <- TRUE
        if (!(identical(intersect(reads1.list, reads2.list), character(0)))) {
            stop(paste("CHIPSTER-NOTE: ", "One or more files is listed in both lists."))
        }
    }
} else if (fileNotOk("reads002.fq")) {
    # Case: No list file, one fastq -> single single-end alignment
    reads1.list <- paste("reads001.fq")
    paired.end <- FALSE
} else if (fileNotOk("reads003.fq")) {
    # Case: No list file, two fastq  -> single paired-end alignment
    in.sorted <- input.names[order(input.names[, 2]), ]
    reads <- grep("reads", in.sorted[, 1], value = TRUE)
    reads1.list <- reads[1]
    reads2.list <- reads[2]
    paired.end <- TRUE
} else {
    # Case: No list files, more than two fastq inputs -> error
    stop(paste("CHIPSTER-NOTE: ", "More than two FASTQ files selected, but no list file(s) provided."))
}

for (i in 1:length(reads1.list)) {
    # command ending
    sam.file <- paste(c(as.character(i), ".sam"), collapse = "")
    bam.file <- paste(c(as.character(i), ".bam"), collapse = "")

    if (paired.end) {
        command.end <- paste(bwa.genome, reads1.list[i], reads2.list[i], "1>", sam.file, "2>> bwa.log")
    } else {
        command.end <- paste(bwa.genome, reads1.list[i], "1>", sam.file, "2>> bwa.log")
    }
    # run bwa alignment
    bwa.command <- paste(command.start, bwa.parameters, command.end)

    documentCommand(bwa.command)

    # stop(paste('CHIPSTER-NOTE: ', bwa.command))
    # why the log is empty if this is run with runExternal()?
    system(bwa.command)

    # convert sam to bam
    runExternal(paste(samtools.binary, "view -b", sam.file, "-o", bam.file))
}

runExternal("echo BWA ready >> bwa.log")
# system("ls -l >> bwa.log")
# Join bam files
if (fileOk("2.bam")) {
    # more than one bam exists, so join them
    runExternal("ls *.bam > bam.list")
    runExternal(paste(samtools.binary, "merge -b bam.list alignment.bam"))
} else {
    # only one bam, so just rename it
    runExternal("mv 1.bam alignment.bam")
}

# Change file named in BAM header to display names
displayNamesToBAM("alignment.bam")

# sort bam
runExternal(paste(samtools.binary, "sort alignment.bam -o alignment.sorted.bam"))

# index bam
runExternal(paste(samtools.binary, "index alignment.sorted.bam"))

# rename result files
runExternal("mv alignment.sorted.bam bwa.bam")
if (index.file == "index_file") {
    runExternal("mv alignment.sorted.bam.bai bwa.bam.bai")
}

# Substitute display names to log for clarity
displayNamesToFile("bwa.log")

# Handle output names
#
# read input names
inputnames <- read_input_definitions()

# Determine base name

# Default name if neither special case match.
# Special cases are for compatibility with older versions of the script
basename <- "bwa_multi"

# Special case 1: Only one input
if (fileNotOk("reads002.fq")) {
    basename <- strip_name(inputnames$reads001.fq)
}
# Special case 2: Paired end and only two inputs
if (paired.end && fileNotOk("reads003.fq")) {
    basename <- paired_name(strip_name(inputnames$reads001.fq), strip_name(inputnames$reads002.fq))
}
# system("ls -l >> bwa.log")
# Make a matrix of output names
outputnames <- matrix(NA, nrow = 3, ncol = 2)
outputnames[1, ] <- c("bwa.bam", paste(basename, ".bam", sep = ""))
outputnames[2, ] <- c("bwa.bam.bai", paste(basename, ".bam.bai", sep = ""))
if (new_index_created == "yes") {
    hg_ifn <- strip_name(inputnames$genome.txt)
    outputnames[3, ] <- c("bwa_index.tar", paste(hg_ifn, "_bwa_index.tar", sep = ""))
}


# Write output definitions file
write_output_definitions(outputnames)

# save version information
bwa.version <- system(paste(bwa.binary, " 2>&1 | grep Version"), intern = TRUE)
documentVersion("BWA", bwa.version)

samtools.version <- system(paste(samtools.binary, "--version | grep samtools"), intern = TRUE)
documentVersion("Samtools", samtools.version)
