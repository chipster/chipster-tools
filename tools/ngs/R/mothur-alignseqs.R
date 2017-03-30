# TOOL mothur-alignseqs.R: "Align sequences to reference" (Given a fasta file of 16S rRNA sequences, this tool aligns them to the Silva reference set available on the server. Alternatively you can give a customized reference fasta. If you do so, make sure the files are correctly assigned in the parameters section! Please note that it can take some time to run this tool. This tool is based on the Mothur tool align.seqs.)
# INPUT reads.fasta: "FASTA file" TYPE FASTA
# INPUT OPTIONAL reference.fasta: "Custom reference FASTA file" TYPE FASTA
# INPUT OPTIONAL a.count_table: "Count table" TYPE MOTHUR_COUNT
# OUTPUT OPTIONAL aligned.fasta.gz
# OUTPUT aligned-summary.tsv
# OUTPUT OPTIONAL custom.reference.summary.tsv
# PARAMETER OPTIONAL reference: "Reference" TYPE [bacterial, full, own] DEFAULT bacterial (Reference sequences to use.)
# PARAMETER OPTIONAL start: "Start" TYPE INTEGER (Start point of your region of interest)
# PARAMETER OPTIONAL end: "End" TYPE INTEGER (End point of your region of interest)
# PARAMETER OPTIONAL keepdots: "Remove leading and trailing dots" TYPE [yes, no] DEFAULT yes (Remove leading and trailing dots.)


# EK 05.06.2013
# ML 21.12.2016 update (new Silva version)
# ML 4.1.2016 new, whole Silva reference
# ML 14.3.2017 reference option (bacterial vs whole)
# ML 15.3.2017 add pcr.seqs options

# OUTPUT log.txt

# check out if the file is compressed and if so unzip it
source(file.path(chipster.common.path, "zip-utils.R"))
unzipIfGZipFile("reads.fasta")

# binary
binary <- c(file.path(chipster.tools.path, "mothur", "mothur"))


if (reference=="bacterial"){
	# new bacterial references:
	data.path <- c(file.path(chipster.tools.path, "mothur-silva-reference", "silva.bacteria"))
	template.path <- c(file.path(data.path, "silva.bacteria.fasta"))
}
if (reference=="full"){
	# new whole references:
	data.path <- c(file.path(chipster.tools.path,"mothur-silva-reference", "mothur-silva-reference-whole"))
	template.path <- c(file.path(data.path, "silva.nr_v123.align")) 
}

# create a symlink, because otherwise the modified reference will go to the reference folder
system(paste("ln -s ", template.path, " template.fasta", sep=""))


# batch file 1 -pcr.seqs -only if user determined end or start 

if (!is.na(start) | !is.na(end)){
	
	pcrseqs.options <- ""
	if (reference=="own"){
		if (file.exists("reference.fasta")) {
		pcrseqs.options <- paste(pcrseqs.options, "pcr.seqs(fasta=reference.fasta", sep="")
		} else{
		stop('CHIPSTER-NOTE: If you choose to use your own reference, you need to give the fasta file for that as input!')
		}
	} else {
		# if using full or bacterial silva reference:
		#pcrseqs.options <- paste(pcrseqs.options, "pcr.seqs(fasta=", template.path, sep="")	
		pcrseqs.options <- paste(pcrseqs.options, "pcr.seqs(fasta=template.fasta", sep="")	
	}

	if (!is.na(start)){
		pcrseqs.options <- paste(pcrseqs.options, ", start=", start, sep="")
	}
	if (!is.na(end)){
		pcrseqs.options <- paste(pcrseqs.options, ", end=", end, sep="")
	}
	if (keepdots=="yes"){
		pcrseqs.options <- paste(pcrseqs.options, ", keepdots=F", sep="")
	}
	if (keepdots=="no"){
		pcrseqs.options <- paste(pcrseqs.options, ", keepdots=T", sep="")
	}
	pcrseqs.options <- paste(pcrseqs.options, ")", sep="")

	# Write batch file
	write(pcrseqs.options, "pcrseq.mth", append=F)
	# command
	command <- paste(binary, "pcrseq.mth", "> log.txt 2>&1")
	# run
	system(command)

}
	
	# rename the ref file as custom.reference.fasta 
	if (file.exists("template.pcr.fasta")) {
			system("mv template.pcr.fasta custom.reference.fasta")
	}	else if (file.exists("reference.pcr.fasta")) {
			system("mv reference.pcr.fasta custom.reference.fasta")
	}	else if (file.exists("template.fasta")) {
			system("mv template.fasta custom.reference.fasta")
	}	
	
	
	
	#  summary file from this step:
	if (file.exists("custom.reference.fasta")) {
		# batch file 2
		write("summary.seqs(fasta=custom.reference.fasta)", "summary.mth", append=F)
		# command
		command2 <- paste(binary, "summary.mth", "> log_raw.txt")
		# run
		system(command2)
		# Post process output
		system("grep -A 10 Start log_raw.txt > custom.reference.summary2.tsv")
		# Remove one tab to get the column naming look nice:
		system("sed 's/^		/	/' custom.reference.summary2.tsv > custom.reference.summary.tsv")
	}	


# batch file 2 -align.seqs
write(paste("align.seqs(fasta=reads.fasta, reference=custom.reference.fasta)", sep=""), "batch.mth", append=F)
# write(paste("align.seqs(fasta=reads.fasta, template=", template.path, ")", sep=""), "batch.mth", append=F)

# command
command <- paste(binary, "batch.mth", ">> log.txt 2>&1")

# run
system(command)

system("mv reads.align aligned.fasta")

# batch file 3 -summary.seqs
if (file.exists("a.count_table")){
	write("summary.seqs(fasta=aligned.fasta, count=a.count_table)", "summary.mth", append=F)
} else {
write("summary.seqs(fasta=aligned.fasta)", "summary.mth", append=F)
}

# command
command2 <- paste(binary, "summary.mth", "> log_raw.txt")

# run
system(command2)

# zip fasta
system("gzip aligned.fasta")

# Post process output
system("grep -A 10 Start log_raw.txt > aligned-summary2.tsv")
# Remove one tab to get the column naming look nice:
system("sed 's/^		/	/' aligned-summary2.tsv > aligned-summary.tsv")
