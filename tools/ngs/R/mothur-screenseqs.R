# TOOL mothur-screenseqs.R: "Screen sequences with Mothur" (Keep sequences that fulfill certain user defined criteria. )
# INPUT reads.trim.unique.fasta: "FASTA file" TYPE FASTA
# INPUT OPTIONAL reads.groups: "Groups file" TYPE MOTHUR_GROUPS
# INPUT OPTIONAL reads.trim.names: "Names file" TYPE MOTHUR_NAMES
# INPUT OPTIONAL reads.summary: "Summary file" TYPE GENERIC
# INPUT OPTIONAL reads.count: "Count file" TYPE GENERIC
# OUTPUT OPTIONAL reads.trim.unique.good.fasta
# OUTPUT OPTIONAL reads.good.groups
# OUTPUT OPTIONAL reads.trim.good.names
# OUTPUT OPTIONAL summary.trim.screen.tsv
# OUTPUT OPTIONAL reads.good.count
# PARAMETER OPTIONAL minlength: "Minimum length of the sequences" TYPE INTEGER (What is the minimum lenght of the sequences?)
# PARAMETER OPTIONAL maxlength: "Maximum length of the sequences" TYPE INTEGER (What is the maximum lenght of the sequences?)
# PARAMETER OPTIONAL end: "End position" TYPE INTEGER (By which position should the sequences end?)
# PARAMETER OPTIONAL start: "Start position" TYPE INTEGER (By which position should the sequences start?)
# PARAMETER OPTIONAL optimize: "Optimize by"  TYPE [empty, minlength, start, end] DEFAULT empty  (Optimize according to minlength, start or end position. Please note that if you use this option, you can't determine the same criteria above! Fill in the optimization criteria below as well.)
# PARAMETER OPTIONAL criteria: "Optimization criteria"  TYPE INTEGER FROM 0 TO 100  (Optimization criteria. For example 85 means that mothur will optimize the cutoff for the above chosen quality so that 85% of the sequences are kept.)
# PARAMETER OPTIONAL maxambig: "Maximum number of ambiguous bases" TYPE INTEGER (How many ambiguous bases are allowed in the sequences?)
# PARAMETER OPTIONAL maxhomop: "Maximum homopolymer length" TYPE INTEGER (Maximum length of homopolymers allowed)


# ML 03.03.2016
#Output File Names: 
#reads.trim.unique.good.fasta
#reads.trim.unique.bad.accnos
#reads.good.groups


# binary
binary <- c(file.path(chipster.tools.path, "mothur", "mothur"))

# Add options
screenseqs.options <- ""
screenseqs.options <- paste(screenseqs.options, "screen.seqs(fasta=reads.trim.unique.fasta")
if (file.exists("reads.trim.names")){
	screenseqs.options <- paste(screenseqs.options, " name=reads.trim.names", sep=",")
}
if (file.exists("reads.groups")){
	screenseqs.options <- paste(screenseqs.options, " group=reads.groups", sep=",")
}
if (file.exists("reads.summary")){
	screenseqs.options <- paste(screenseqs.options, " summary=reads.summary", sep=",")
}
if (file.exists("reads.count")){
	screenseqs.options <- paste(screenseqs.options, " count=reads.count", sep=",")
}
# Sanity check (User can't optimize by minlength and specify a minlength at the same time)
if (optimize != "empty"){
	if ( (optimize == "minlength" && !is.na(minlength)) || (optimize == "start" && !is.na(start)) || (optimize == "end" && !is.na(end))){
		stop('CHIPSTER-NOTE: You cant determine minlenght and choose to optimize according the same criteria: choose only one of these!')
	}
}	
	
if (!is.na(minlength)){
	screenseqs.options <- paste(screenseqs.options, ", minlength=", minlength, sep="")
}
if (!is.na(maxlength)){
	screenseqs.options <- paste(screenseqs.options, ", maxlength=", maxlength, sep="")
}
if (!is.na(end)){
	screenseqs.options <- paste(screenseqs.options, ", end=", end, sep="")
}
if (!is.na(start)){
	screenseqs.options <- paste(screenseqs.options, ", start=", start, sep="")
}
if (!is.na(maxambig)){
	screenseqs.options <- paste(screenseqs.options, ", maxambig=", maxambig, sep="")
}
if (!is.na(maxhomop)){
	screenseqs.options <- paste(screenseqs.options, ", maxhomop=", maxhomop, sep="")
}
if (optimize != "empty"){
	screenseqs.options <- paste(screenseqs.options, ", optimize=", optimize, sep="")
	screenseqs.options <- paste(screenseqs.options, ", criteria=", criteria, sep="")
}

screenseqs.options <- paste(screenseqs.options, ")", sep="")

# stop(paste('CHIPSTER-NOTE: ', screenseqs.options))

# Write batch file
write(screenseqs.options, "trim.mth", append=F)
#write("screen.seqs(fasta=reads.trim.fasta)", "trim.mth", append=T)

# command
command <- paste(binary, "trim.mth", "> log.txt 2>&1")


# run
system(command)

# rename the file
system("mv reads.trim.unique.good.align reads.trim.unique.good.fasta")

# batch file
# write("summary.seqs(fasta=reads.trim.unique.good.fasta, name=reads.trim.good.names)", "summary.mth", append=F)
write("summary.seqs(fasta=reads.trim.unique.good.fasta)", "summary.mth", append=F)

# command
command <- paste(binary, "summary.mth", "> log_raw.txt")

# run
system(command)

# Postprocess output files
system("grep -A 10 Start log_raw.txt > summary.trim.screen2.tsv")
# Remove one tab to get the column naming look nice:
system("sed 's/^		/	/' summary.trim.screen2.tsv > summary.trim.screen.tsv")

