# TOOL mothur-screenseqs.R: "Screen sequences with Mothur" (Keep sequences that fulfill certain user defined criteria. )
# INPUT reads.trim.unique.fasta: "FASTA file" TYPE FASTA
# INPUT reads.groups: "Groups file" TYPE MOTHUR_GROUPS
# INPUT reads.trim.names: "Names file" TYPE MOTHUR_NAMES
# OUTPUT OPTIONAL reads.trim.unique.good.fasta
# OUTPUT OPTIONAL reads.good.groups
# OUTPUT OPTIONAL reads.trim.good.names
# OUTPUT OPTIONAL summary.trim.screen.tsv
# PARAMETER OPTIONAL minlength: "Minimum length of the sequences" TYPE INTEGER (How long should the sequences at least be?)
# PARAMETER OPTIONAL end: "End position" TYPE INTEGER (By which position should the sequences end?)
# PARAMETER OPTIONAL start: "Start position" TYPE INTEGER (By which position should the sequences start?)
# PARAMETER OPTIONAL optimize: "Optimize by"  TYPE [empty, minlength, start, end] DEFAULT empty  (Optimize according to minlength, start or end position. Please note that if you use this option, you can't determine the same criteria above! Fill in the optimization criteria below as well.)
# PARAMETER OPTIONAL criteria: "Optimization criteria"  TYPE INTEGER FROM 0 TO 100  (Optimization criteria. For example 85 means that mothur will optimize the cutoff for the above chosen quality so that 85% of the sequences are kept.)



# ML 03.03.2016

# binary
binary <- c(file.path(chipster.tools.path, "mothur", "mothur"))

# Add options
screenseqs.options <- ""
screenseqs.options <- paste(screenseqs.options, "screen.seqs(fasta=reads.trim.unique.fasta, name=reads.trim.names, group=reads.groups")

# Sanity check (User can't optimize by minlength and specify a minlength at the same time)
if (optimize != "empty"){
	if ( (optimize == "minlength" && !is.na(minlength)) || (optimize == "start" && !is.na(start)) || (optimize == "end" && !is.na(end))){
		stop('CHIPSTER-NOTE: You cant determine minlenght and choose to optimize according the same criteria: choose only one of these!')
	}
}	
	
if (!is.na(minlength)){
	screenseqs.options <- paste(screenseqs.options, ", minlength=", minlength, sep="")
}
if (!is.na(end)){
	screenseqs.options <- paste(screenseqs.options, ", end=", end, sep="")
}
if (!is.na(start)){
	screenseqs.options <- paste(screenseqs.options, ", start=", start, sep="")
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
command <- paste(binary, "trim.mth")

# run
system(command)

# rename the file
system("mv reads.trim.unique.good.align reads.trim.unique.good.fasta")

# batch file
write("summary.seqs(fasta=reads.trim.unique.good.fasta, name=reads.trim.good.names)", "summary.mth", append=F)

# command
command <- paste(binary, "summary.mth", "> log_raw.txt")

# run
system(command)

# Postprocess output files
system("grep -A 9 Start log_raw.txt > summary.trim.screen2.tsv")
# Remove one tab to get the column naming look nice:
system("sed 's/^		/	/' summary.trim.screen2.tsv > summary.trim.screen.tsv")

