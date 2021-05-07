# TOOL mothur-screenseqs.R: "Screen sequences for several criteria" (Keeps sequences that fulfill user-defined criteria. This tool is based on the Mothur tool screen.seqs.)
# INPUT a.fasta: "FASTA file" TYPE GENERIC
# INPUT OPTIONAL a.groups: "Groups file" TYPE MOTHUR_GROUPS
# INPUT OPTIONAL a.count: "Count file" TYPE MOTHUR_COUNT
# OUTPUT OPTIONAL screened.fasta.gz
# OUTPUT OPTIONAL screened.groups
# OUTPUT OPTIONAL summary.screened.tsv
# OUTPUT OPTIONAL screened.count_table
# PARAMETER OPTIONAL maxambig: "Maximum number of ambiguous bases" TYPE INTEGER (How many ambiguous bases are allowed in a sequence)
# PARAMETER OPTIONAL maxhomop: "Maximum homopolymer length" TYPE INTEGER (Maximum length of homopolymers allowed)
# PARAMETER OPTIONAL minlength: "Minimum length" TYPE INTEGER (What is the minimum length of the sequences to be kept?)
# PARAMETER OPTIONAL maxlength: "Maximum length" TYPE INTEGER (What is the maximum length of the sequences to be kept?)
# PARAMETER OPTIONAL start: "Alignment start position" TYPE INTEGER (Remove sequences which start after this position)
# PARAMETER OPTIONAL end: "Alignment end position" TYPE INTEGER (Remove sequences which end before this position)
# PARAMETER OPTIONAL optimize: "Optimize by"  TYPE [empty, minlength, start, end] DEFAULT empty  (Optimize according to minlength, start or end position. Please note that if you use this option, you can't determine the same criteria above! Fill in the optimization criteria below as well.)
# PARAMETER OPTIONAL criteria: "Optimization criteria"  TYPE INTEGER FROM 0 TO 100  (Optimization criteria. For example 85 means that Mothur will optimize the cutoff for the above chosen quality so that 85% of the sequences are kept.)


# ML 03.03.2016
# ML 17.3.2017 Clarify inputs and outputs
#Output File Names: 
#reads.trim.unique.good.fasta
#reads.trim.unique.bad.accnos
#reads.good.groups
# OUTPUT OPTIONAL screened.names
# INPUT OPTIONAL a.names: "Names file" TYPE MOTHUR_NAMES

source(file.path(chipster.common.path,"tool-utils.R"))
source(file.path(chipster.common.path,"zip-utils.R"))

# Check if fasta is zipped and unzip it if needed
unzipIfGZipFile("a.fasta")

# binary
binary <- c(file.path(chipster.tools.path,"mothur","mothur"))
version <- system(paste(binary,"--version"),intern = TRUE)
documentVersion("Mothur",version)

# Add options
screenseqs.options <- ""
screenseqs.options <- paste(screenseqs.options,"screen.seqs(fasta=a.fasta")
if (file.exists("a.names")) {
  screenseqs.options <- paste(screenseqs.options," name=a.names",sep = ",")
}
if (file.exists("a.groups")) {
  screenseqs.options <- paste(screenseqs.options," group=a.groups",sep = ",")
}
if (file.exists("a.summary")) {
  screenseqs.options <- paste(screenseqs.options," summary=a.summary",sep = ",")
}
if (file.exists("a.count")) {
  screenseqs.options <- paste(screenseqs.options," count=a.count",sep = ",")
}
# Sanity check (User can't optimize by minlength and specify a minlength at the same time)
if (optimize != "empty") {
  if ((optimize == "minlength" && !is.na(minlength)) || (optimize == "start" && !is.na(start)) || (optimize == "end" && !is.na(end))) {
    stop('CHIPSTER-NOTE: You cant determine minlenght and choose to optimize according the same criteria: choose only one of these!')
  }
}

if (!is.na(minlength)) {
  screenseqs.options <- paste(screenseqs.options,", minlength=",minlength,sep = "")
}
if (!is.na(maxlength)) {
  screenseqs.options <- paste(screenseqs.options,", maxlength=",maxlength,sep = "")
}
if (!is.na(end)) {
  screenseqs.options <- paste(screenseqs.options,", end=",end,sep = "")
}
if (!is.na(start)) {
  screenseqs.options <- paste(screenseqs.options,", start=",start,sep = "")
}
if (!is.na(maxambig)) {
  screenseqs.options <- paste(screenseqs.options,", maxambig=",maxambig,sep = "")
}
if (!is.na(maxhomop)) {
  screenseqs.options <- paste(screenseqs.options,", maxhomop=",maxhomop,sep = "")
}
if (optimize != "empty") {
  screenseqs.options <- paste(screenseqs.options,", optimize=",optimize,sep = "")
  screenseqs.options <- paste(screenseqs.options,", criteria=",criteria,sep = "")
}
screenseqs.options <- paste(screenseqs.options,", processors=",chipster.threads.max,sep = "")
screenseqs.options <- paste(screenseqs.options,")",sep = "")

# stop(paste('CHIPSTER-NOTE: ', screenseqs.options))

# Write batch file
documentCommand(screenseqs.options)
write(screenseqs.options,"trim.mth",append = FALSE)
#write("screen.seqs(fasta=reads.trim.fasta)", "trim.mth", append=T)

# command
command <- paste(binary,"trim.mth","> log.txt 2>&1")

# run
system(command)

# rename the result files
system("mv a.good.fasta screened.fasta")
if (file.exists("a.good.count")) {
  system("mv a.good.count screened.count_table")
}
if (file.exists("a.good.groups")) {
  system("mv a.good.groups screened.groups")
}

# batch file
# write("summary.seqs(fasta=reads.trim.unique.good.fasta, name=reads.trim.good.names)", "summary.mth", append=F)
# write("summary.seqs(fasta=screened.fasta)", "summary.mth", append=F)
summaryseqs.options <- paste("summary.seqs(fasta=screened.fasta")
if (file.exists("screened.count_table")) {
  summaryseqs.options <- paste(summaryseqs.options,", count=screened.count_table")
}
summaryseqs.options <- paste(summaryseqs.options,", processors=",chipster.threads.max,")",sep = "")
documentCommand(summaryseqs.options)
write(summaryseqs.options,"summary.mth",append = FALSE)
# command
command <- paste(binary,"summary.mth","> log_raw.txt")

# run
system(command)

# zip output fasta
system("gzip screened.fasta")

# Postprocess output files
system("grep -A 10 Start log_raw.txt > summary.screen2.tsv")
# Remove one tab to get the column naming look nice:
system("sed 's/^		/	/' summary.screen2.tsv > summary.screened.tsv")

