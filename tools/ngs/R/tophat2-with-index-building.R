# TOOL tophat2-with-index-building.R: "TopHat2 for paired end reads and own genome" (This tool uses TopHat2 to align paired-end reads to a provided reference genome to identify exon-exon splice junctions. You need to supply the paired-end reads in FASTQ format. You also need to provide the TopHat2 index for your reference genome in TAR format. )
# INPUT reads{...}.fq: "Reads" TYPE FASTQ
# INPUT OPTIONAL reads1.txt: "List of read 1 files" TYPE TEXT_STRICT
# INPUT OPTIONAL reads2.txt: "List of read 2 files" TYPE TEXT_STRICT
# INPUT OPTIONAL genome.tar: "Tophat2 index of the genome to align against" TYPE TAR
# OUTPUT OPTIONAL tophat.bam
# OUTPUT OPTIONAL tophat.bam.bai
# OUTPUT OPTIONAL junctions.bed
# OUTPUT OPTIONAL tophat-summary.txt
# OUTPUT OPTIONAL logs.tar
# PARAMETER library.type: "Library type" TYPE [fr-unstranded: fr-unstranded, fr-firststrand: fr-firststrand, fr-secondstrand: fr-secondstrand] DEFAULT fr-unstranded (Which library type to use. For directional\/strand specific library prepartion methods, choose fr-firststrand or fr-secondstrand depending on the preparation method: if the first read \(read1\) maps to the opposite, non-coding strand, choose fr-firststrand. If the first read maps to the coding strand, choose fr-secondstrand. For example for Illumina TruSeq Stranded sample prep, choose fr-firstsrand.)
# PARAMETER  mate.inner.distance: "Expected inner distance between mate pairs" TYPE INTEGER DEFAULT 200 (Expected mean inner distance between mate pairs. For example, if your fragment size is 300 bp and read length is 50 bp, the inner distance is 200.)
# PARAMETER OPTIONAL no.novel.juncs: "When transcriptome index is included in the index package, ignore novel junctions" TYPE [yes, no] DEFAULT no (Only look for reads across junctions indicated in the supplied transcriptome index.)
# PARAMETER OPTIONAL quality.format: "Base quality encoding used" TYPE [sanger: "Sanger - Phred+33", phred64: "Phred+64"] DEFAULT sanger (Quality encoding used in the fastq file.)
# PARAMETER OPTIONAL mate.std.dev: "Standard deviation for the inner distances between mate pairs" TYPE INTEGER DEFAULT 20 (The standard deviation for the distribution on inner distances between mate pairs. The default is 20bp.)
# PARAMETER OPTIONAL max.multihits: "How many hits is a read allowed to have" TYPE INTEGER FROM 1 TO 1000000 DEFAULT 20 (Instructs TopHat to allow up to this many alignments to the reference for a given read.)
# PARAMETER OPTIONAL mismatches: "Number of mismatches allowed in final alignment" TYPE INTEGER FROM 0 TO 5 DEFAULT 2 (Final read alignments having more than this many mismatches are discarded.)
# PARAMETER OPTIONAL min.anchor.length: "Minimum anchor length" TYPE INTEGER FROM 3 TO 1000 DEFAULT 8 (TopHat2 will report junctions spanned by reads with at least this many bases on each side of the junction. Note that individual spliced alignments may span a junction with fewer than this many bases on one side. However, every junction involved in spliced alignments is supported by at least one read with this many bases on each side.)
# PARAMETER OPTIONAL splice.mismatches: "Maximum number of mismatches allowed in the anchor" TYPE INTEGER FROM 0 TO 2 DEFAULT 0 (The maximum number of mismatches that may appear in the anchor region of a spliced alignment.)
# PARAMETER OPTIONAL min.intron.length: "Minimum intron length" TYPE INTEGER FROM 4 TO 1000 DEFAULT 70 (TopHat2 will ignore donor-acceptor pairs closer than this many bases apart.)
# PARAMETER OPTIONAL max.intron.length: "Maximum intron length" TYPE INTEGER FROM 1 TO 1000000 DEFAULT 500000 (TopHat2 will ignore donor-acceptor pairs farther than this many bases apart, except when such a pair is supported by a split segment alignment of a long read.)
# PARAMETER OPTIONAL no.mixed: "Report only paired alignments" TYPE [yes, no] DEFAULT yes (Only report read alignments if both reads in a pair can be mapped.)


# EK 17.4.2012 added -G and -g options
# MG 24.4.2012 added ability to use gtf files from Chipster server
# AMS 19.6.2012 added unzipping
# AMS 27.6.2012 added parameter mate.std.dev, allow negative values for mate.inner.distance
# AMS 4.10.2012 added BED sorting
# KM 10.7. 2012 added RN5
# AMS 11.11.2013 added thread support
# AMS 3.1.2014 added transcriptome index for human
# EK 3.1.2014 added alignment summary to output, added quality and mismatch parameter
# AMS 22.5.2014 modified to use own genome
# ML 15.01.2015 Added the library-type parameter
# AMS 29.01.2015 Removed optional outputs deletions.bed and insertions.bed

# PARAMETER OPTIONAL no.discordant: "Report only concordant alignments" TYPE [yes, no] DEFAULT yes (Report only concordant mappings.)

source(file.path(chipster.common.lib.path, "tool-utils.R"))
source(file.path(chipster.common.lib.path, "zip-utils.R"))

# check out if the file is compressed and if so unzip it
input.names <- read.table("chipster-inputs.tsv", header = F, sep = "\t")
for (i in 1:nrow(input.names)) {
    unzipIfGZipFile(input.names[i, 1])
}



options(scipen = 10)
# max.intron.length <- formatC(max.intron.length, "f", digits = 0)

# setting up TopHat
tophat.binary <- c(file.path(chipster.tools.path, "tophat2", "tophat2"))
version <- system(paste(tophat.binary, "--version | cut -d ' ' -f 2"), intern = TRUE)
documentVersion("TopHat", version)

bowtie.binary <- c(file.path(chipster.tools.path, "bowtie2-2.2.9", "bowtie2"))
version <- system(paste(bowtie.binary, "--version | head -1 | cut -d ' ' -f 3"), intern = TRUE)
documentVersion("Bowtie", version)
bowtie2.index.binary <- file.path(chipster.module.path, "shell", "check_bowtie2_index.sh")

path.bowtie <- c(file.path(chipster.tools.path, "bowtie2"))
path.samtools <- c(file.path(chipster.tools.path, "samtools-0.1.19"))
set.path <- paste(sep = "", "PATH=", path.bowtie, ":", path.samtools, ":$PATH")

print("Extracting tar formatted gemome index file")
runExternal("tar xf genome.tar")

# Check bowtie2 index base name
if (length(Sys.glob("*.1.bt2")) != 0) {
    f <- list.files(getwd(), pattern = "\\.1.bt2$")
    bowtie2.genome <- substr(f[1], 1, nchar(f[1]) - 6)
} else {
    stop("CHIPSTER-NOTE: The .tar package does not seem to contain a valid TopHat2 index.")
}

# Check optional tophat2 index base name
if (length(Sys.glob("tophat2/*.1.bt2")) != 0) {
    f <- list.files(file.path(getwd(), "tophat2"), pattern = "\\.1.bt2$")
    tophat2.genome <- file.path("tophat2", substr(f[1], 1, nchar(f[1]) - 6))
    print(paste("TopHat2 genome name:", tophat2.genome))
} else {
    print("The .tar package does not contain TopHat2 transcriptome index.")
    tophat2.genome <- NULL
}


runExternal("ls -lah")

print(paste("Bowtie2 genome name:", bowtie2.genome))

# command start
command.start <- paste("bash -c '", set.path, tophat.binary)

# parameters
# command.parameters <- paste("--bowtie1 -r", mate.inner.distance, "--mate-std-dev", mate.std.dev, "-a", min.anchor.length, "-m", splice.mismatches, "-i", min.intron.length, "-I", max.intron.length, "-g", max.multihits, "--library-type fr-unstranded")
# command.parameters <- paste("-p", chipster.threads.max, "-r", mate.inner.distance, "--mate-std-dev", mate.std.dev, "--read-mismatches", mismatches, "-a", min.anchor.length, "-m", splice.mismatches, "-i", min.intron.length, "-I", max.intron.length, "-g", max.multihits, "--library-type fr-unstranded")
command.parameters <- paste("-p", chipster.threads.max, "-r", mate.inner.distance, "--mate-std-dev", mate.std.dev, "--read-mismatches", mismatches, "-a", min.anchor.length, "-m", splice.mismatches, "-i", min.intron.length, "-I", max.intron.length, "-g", max.multihits, "--library-type", library.type)

if (mismatches > 2) {
    command.parameters <- paste(command.parameters, "--read-edit-dist", mismatches)
}

if (quality.format == "phred64") {
    command.parameters <- paste(command.parameters, "--phred64-quals")
}

# if (no.discordant == "yes"){
# 	command.parameters <- paste(command.parameters, "--no-discordant")
# }

if (no.mixed == "yes") {
    command.parameters <- paste(command.parameters, "--no-mixed")
}

# optional transriptome index
if (!is.null(tophat2.genome)) {
    command.parameters <- paste(command.parameters, "--transcriptome-index", tophat2.genome)
    if (no.novel.juncs == "yes") {
        command.parameters <- paste(command.parameters, "--no-novel-juncs")
    }
}


# Input files
if (file.exists("reads1.txt") && file.exists("reads2.txt")) {
    # Case: list files exist
    reads1.list <- make_input_list("reads1.txt")
    reads2.list <- make_input_list("reads2.txt")
    if (identical(intersect(reads1.list, reads2.list), character(0))) {
        reads1 <- paste(reads1.list, sep = "", collapse = ",")
        reads2 <- paste(reads2.list, sep = "", collapse = ",")
    } else {
        stop(paste("CHIPSTER-NOTE: ", "One or more files is listed in both lists."))
    }
} else if (file.exists("reads002.fq") && !file.exists("reads003.fq")) {
    # Case: no list file, but only two fastq inputs
    in.sorted <- input.names[order(input.names[, 2]), ]
    reads <- grep("reads", in.sorted[, 1], value = TRUE)
    reads1 <- reads[1]
    reads2 <- reads[2]
} else {
    # Case: no list files, more than two fastq inputs
    stop(paste("CHIPSTER-NOTE: ", "List file is missing. You need to provide a list of read files for both directions."))
}

# command ending
# tophat2 prints status to stderr
command.end <- paste(bowtie2.genome, reads1, reads2, "2>&1'")

# run tophat
command <- paste(command.start, command.parameters, command.end)

documentCommand(command)

runExternal(command)


# samtools binary
samtools.binary <- c(file.path(chipster.tools.path, "samtools-0.1.19", "samtools"))
version <- system(paste(samtools.binary, "2>&1 |head -3 |tail -1 | cut -d ' ' -f 2"), intern = TRUE)
documentVersion("SAMtools", version)

# sort bam (removed because TopHat itself does the sorting)
# system(paste(samtools.binary, "sort tophat_out/accepted_hits.bam tophat"))
runExternal("mv tophat_out/accepted_hits.bam tophat.bam")

# index bam
runExternal(paste(samtools.binary, "index tophat.bam"))

runExternal("mv tophat_out/junctions.bed junctions.u.bed")
runExternal("mv tophat_out/insertions.bed insertions.u.bed")
runExternal("mv tophat_out/deletions.bed deletions.u.bed")
runExternal("mv tophat_out/align_summary.txt tophat-summary.txt")

# sorting BEDs
source(file.path(chipster.common.lib.path, "bed-utils.R"))

no.results <- "TRUE"

if (file.exists("junctions.u.bed")) {
    size <- file.info("junctions.u.bed")$size
    if (size > 100) {
        bed <- read.table(file = "junctions.u.bed", skip = 1, sep = "\t")
        colnames(bed)[1:2] <- c("chr", "start")
        sorted.bed <- sort.bed(bed)
        write.table(sorted.bed, file = "junctions.bed", sep = "\t", row.names = F, col.names = F, quote = F)
        no.results <- "FALSE"
    }
}

if (file.exists("insertions.u.bed")) {
    size <- file.info("insertions.u.bed")$size
    if (size > 100) {
        bed <- read.table(file = "insertions.u.bed", skip = 1, sep = "\t")
        colnames(bed)[1:2] <- c("chr", "start")
        sorted.bed <- sort.bed(bed)
        write.table(sorted.bed, file = "insertions.bed", sep = "\t", row.names = F, col.names = F, quote = F)
        no.results <- "FALSE"
    }
}

if (file.exists("deletions.u.bed")) {
    size <- file.info("deletions.u.bed")$size
    if (size > 100) {
        bed <- read.table(file = "deletions.u.bed", skip = 1, sep = "\t")
        colnames(bed)[1:2] <- c("chr", "start")
        sorted.bed <- sort.bed(bed)
        write.table(sorted.bed, file = "deletions.bed", sep = "\t", row.names = F, col.names = F, quote = F)
        no.results <- "FALSE"
    }
}

# If no BAM file is produced, return the whole logs folder as a tar package
if (fileNotOk("tophat.bam")) {
    runExternal("tar cf logs.tar tophat_out/logs/*")
}

# Handle output names
source(file.path(chipster.common.lib.path, "tool-utils.R"))

# read input names
inputnames <- read_input_definitions()

# Determine base name
name1 <- unlist(strsplit(reads1, ","))
base1 <- strip_name(inputnames[[name1[1]]])

name2 <- unlist(strsplit(reads2, ","))
base2 <- strip_name(inputnames[[name2[1]]])

basename <- paired_name(base1, base2)

# Make a matrix of output names
outputnames <- matrix(NA, nrow = 2, ncol = 2)
outputnames[1, ] <- c("tophat.bam", paste(basename, ".bam", sep = ""))
outputnames[2, ] <- c("tophat.bam.bai", paste(basename, ".bam.bai", sep = ""))

# Write output definitions file
write_output_definitions(outputnames)

# EOF
