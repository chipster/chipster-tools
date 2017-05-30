# TOOL mothur-chimera-uchime.R: "Remove chimeric sequences" (Remove chimeric sequences from a fasta-formatted alignment using the uchime method and the 16S rRNA Silva gold reference set. This tool is based on the Mothur package. Please note that it can take some time to run this tool!)
# INPUT a.fasta: "FASTA file" TYPE FASTA
# INPUT OPTIONAL a.groups: "Groups file" TYPE MOTHUR_GROUPS
# INPUT OPTIONAL a.names: "Names file" TYPE MOTHUR_NAMES
# INPUT OPTIONAL a.count_table: "Count table" TYPE GENERIC
# OUTPUT OPTIONAL chimeras.removed.fasta
# OUTPUT OPTIONAL chimeras.removed.summary.tsv
# OUTPUT OPTIONAL chimeras.removed.count_table
# PARAMETER OPTIONAL dereplicate: "Dereplicate" TYPE [false, true] DEFAULT false (False = If sequence is flagged as chimeric, remove it from all samples)
# PARAMETER OPTIONAL reference: "Reference" TYPE [bacterial, full, none] DEFAULT bacterial (Reference sequences to use. If you choose none, note that you have to give count table as input.)


# OUTPUT OPTIONAL log.txt
# OUTPUT OPTIONAL log2.txt

# EK 18.06.2013
# 30.3.2016 ML added input options and parameters
# ML 21.12.2016 update (new version, new Silva version)
# ML 14.3.2017 reference option (bacterial vs whole) + count-file output

# check out if the file is compressed and if so unzip it
source(file.path(chipster.common.path, "zip-utils.R"))
unzipIfGZipFile("a.fasta")

# binary
binary <- c(file.path(chipster.tools.path, "mothur", "mothur"))
# old references:
#data.path <- c(file.path(chipster.tools.path, "mothur-data"))
#template.path <- c(file.path(data.path, "silva.gold.align"))

uchime.options <- ""

if (reference=="bacterial"){
	# new bacterial references:
	data.path <- c(file.path(chipster.tools.path, "mothur-silva-reference", "silva.bacteria"))
	template.path <- c(file.path(data.path, "silva.gold.ng.fasta"))
	uchime.options <- paste(uchime.options, "chimera.uchime(fasta=a.fasta, template=", template.path,sep="")
}
if (reference=="full"){
	# new whole references:
	data.path <- c(file.path(chipster.tools.path,"mothur-silva-reference", "mothur-silva-reference-whole"))
	template.path <- c(file.path(data.path, "silva.gold.align")) 
	uchime.options <- paste(uchime.options, "chimera.uchime(fasta=a.fasta, template=", template.path,sep="")
}
if (reference=="none"){
	uchime.options <- paste(uchime.options, "chimera.uchime(fasta=a.fasta", sep="")
	if (file.exists("a.count_table")){
		uchime.options <- paste(uchime.options, " count=a.count_table", sep=",")
	}
	else {
		stop('CHIPSTER-NOTE: In order to use the denovo option, note that you have to give a count table as input!')
	}
}
	

if (file.exists("a.names")){
	uchime.options <- paste(uchime.options, " name=a.names", sep=",")
}
if (file.exists("a.groups")){
	uchime.options <- paste(uchime.options, " group=a.groups", sep=",")
}

if (dereplicate=="false"){
	uchime.options <- paste(uchime.options, ", dereplicate=F", sep="")
}
if (dereplicate=="true"){
	uchime.options <- paste(uchime.options, ", dereplicate=T", sep="")
}

uchime.options <- paste(uchime.options, ", processors=", chipster.threads.max, sep="")

uchime.options <- paste(uchime.options, ")", sep="")

#stop(paste('CHIPSTER-NOTE: ', uchime.options))

# Write batch file
write(uchime.options, "batch.mth", append=F)
#write(paste("chimera.uchime(fasta=a.fasta, template=", template.path, ")", sep=""), "batch.mth", append=F)

# command
command <- paste(binary, "batch.mth", "> log.txt 2>&1")

# run
system(command)
#stop(paste('CHIPSTER-NOTE: ', uchime.options))
#Output File Names: 
#a.ref.uchime.chimeras
#a.ref.uchime.accnos
#a.denovo.uchime.chimeras
#a.denovo.uchime.accnos

# batch file 2

if (reference=="none"){
	write("remove.seqs(accnos=a.denovo.uchime.accnos, fasta=a.fasta, count=a.count_table)", "remove.mth", append=F)
} else {
	if(file.exists("a.count_table")){
		write("remove.seqs(accnos=a.ref.uchime.accnos, fasta=a.fasta, count=a.count_table)", "remove.mth", append=F)
	}else {
		write("remove.seqs(accnos=a.ref.uchime.accnos, fasta=a.fasta)", "remove.mth", append=F)
	}
}

# command
command2 <- paste(binary, "remove.mth", "> log2.txt 2>&1")

# run
system(command2)

# Post process output
system("mv a.pick.fasta chimeras.removed.fasta")
if (file.exists("a.pick.count_table")){
	system("mv a.pick.count_table chimeras.removed.count_table") 
	write("summary.seqs(fasta=chimeras.removed.fasta, count=chimeras.removed.count_table)", "summary.mth", append=F)
} else {
	write("summary.seqs(fasta=chimeras.removed.fasta)", "summary.mth", append=F)
}

# command 3
command3 <- paste(binary, "summary.mth", "> log_raw.txt")

# run
system(command3)

# Post process output
system("grep -A 10 Start log_raw.txt > chimeras.removed.summary2.tsv")
# Remove one tab to get the column naming look nice:
system("sed 's/^		/	/' chimeras.removed.summary2.tsv > chimeras.removed.summary.tsv")