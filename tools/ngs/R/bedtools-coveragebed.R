# TOOL bedtools-coveragebed.R: "Coverage BED" (Returns the depth and breadth of coverage of features from A on the intervals in B. This tool is based on the BEDTools package.)
# INPUT file.a: "Input file A" TYPE GENERIC
# INPUT file.b: "Input file B" TYPE GENERIC
# OUTPUT OPTIONAL coveragebed.bed
# OUTPUT OPTIONAL error.txt
# PARAMETER abam: "File A is BAM format" TYPE [yes, no] DEFAULT no (Select yes if file A is BAM format.)
# PARAMETER OPTIONAL s: "Force strandedness" TYPE [yes, no] DEFAULT no (Force strandedness. That is, only include hits in A that overlap B on the same strand. By default, hits are included without respect to strand.)
# PARAMETER OPTIONAL hist: "Report a histogram of coverage" TYPE [yes, no] DEFAULT no (Report a histogram of coverage for each feature in B as well as a summary histogram for all features in B. Output (tab delimited\) after each feature in B: 1\) depth, 2\) \# bases at depth, 3\) size of B, 4\) % of B at depth.)
# PARAMETER OPTIONAL d: "Report the depth at each position" TYPE [yes, no] DEFAULT no (Report the depth at each position in each B feature. Positions reported are one based.  Each position and depth follow the complete B feature.)
# PARAMETER OPTIONAL split: "Treat split BAM or BED12 entries as distinct BED intervals" TYPE [yes, no] DEFAULT no (Treat "split" BAM or BED12 entries as distinct BED intervals when computing coverage. For BAM files, this uses the CIGAR N and D operations to infer the blocks for computing coverage. For BED12 files, this uses the BlockCount, BlockStarts, and BlockEnds fields (i.e., columns 10,11,12\).)

# AMS 23.4.2012
# AMS 11.10.2012 Fixed BAM file support
# AMS 23.9.2013 Improved output/error file handling

# binary
binary <- c(file.path(chipster.tools.path, "bedtools", "bin", "coverageBed"))

# optional options
options <- paste("")
if (s == "yes") {
    options <- paste(options, "-s")
}
if (hist == "yes") {
    options <- paste(options, "-hist")
}
if (d == "yes") {
    options <- paste(options, "-d")
}
if (split == "yes") {
    options <- paste(options, "-split")
}

# input files
if (abam == "yes") {
    options <- paste(options, "-abam file.a -b file.b")
}
if (abam == "no") {
    options <- paste(options, "-a file.a -b file.b")
}

# command
command <- paste(binary, options, "> coveragebed.tmp 2> error.tmp")

# run
system(command)

# Generate output/error message
if (file.info("coveragebed.tmp")$size > 0) {
    system("mv coveragebed.tmp coveragebed.bed")
} else if (file.info("error.tmp")$size > 0) {
    system("mv error.tmp error.txt")
} else {
    system("echo \"# No results found\" > error.txt")
}
