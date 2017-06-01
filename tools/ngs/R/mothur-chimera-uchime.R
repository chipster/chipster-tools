# TOOL mothur-chimera-uchime.R: "Remove chimeric sequences" (Removes chimeric sequences from a fasta-formatted alignment using the vsearch or uchime method. You can select de novo or reference based chimera detection. As a reference you can use the bacterial subset of the 16S rRNA Silva gold reference set. For the uchime method also the full 16S rRNA Silva gold reference set is available. This tool is based on the Mothur tools Chimera.vsearch, Chimera.uchime and Remove.seqs. Please note that it can take some time to run this tool!)
# INPUT a.fasta: "FASTA file" TYPE FASTA
# INPUT OPTIONAL a.count_table: "Count table" TYPE GENERIC
# OUTPUT OPTIONAL chimeras.removed.fasta
# OUTPUT OPTIONAL chimeras.removed.summary.tsv
# OUTPUT OPTIONAL chimeras.removed.count_table
# PARAMETER OPTIONAL reference: "Reference" TYPE [bacterial, full, none] DEFAULT bacterial (Reference sequences to use. Note that if you choose none, you have to give count table as input.)
# PARAMETER OPTIONAL method: "Method" TYPE [vsearch, uchime] DEFAULT vsearch (Chimera detection method to use. Note that vsearch is much faster than the Uchime method.)
# PARAMETER OPTIONAL dereplicate: "Dereplicate" TYPE [false, true] DEFAULT false (False = If a sequence is flagged as chimeric in one sample, remove it from all samples.)


# OUTPUT OPTIONAL log.txt
# OUTPUT OPTIONAL log2.txt

# EK 18.06.2013
# 30.3.2016 ML added input options and parameters
# ML 21.12.2016 update (new version, new Silva version)
# ML 14.3.2017 reference option (bacterial vs whole) + count-file output
# AMS 30.5.2017 added the possibility to use more processors
# EK 1.6.2017 added the vsearch method. Removed the input option for names and groups file as we use the more compact count file now for duplicates.

# INPUT OPTIONAL a.names: "Names file" TYPE MOTHUR_NAMES
# INPUT OPTIONAL a.groups: "Groups file" TYPE MOTHUR_GROUPS

# check out if the file is compressed and if so unzip it
source(file.path(chipster.common.path, "zip-utils.R"))
unzipIfGZipFile("a.fasta")

# binary
binary <- c(file.path(chipster.tools.path, "mothur", "mothur"))
# old references:
#data.path <- c(file.path(chipster.tools.path, "mothur-data"))
#template.path <- c(file.path(data.path, "silva.gold.align"))

# start building the command for the Mothur uchime and vsearch tools based on the method chosen. Set also the reference parameter name based on the method.
if (method=="vsearch"){
	chimera.options <- paste("chimera.vsearch(fasta=a.fasta, processors=", chipster.threads.max, ",", sep="")
	ref <- "reference"
}
if (method=="uchime"){
	chimera.options <- paste("chimera.uchime(fasta=a.fasta, processors=", chipster.threads.max, ",", sep="")
	ref <- "template"
}

if (reference=="bacterial"){
	# bacterial reference
	data.path <- c(file.path(chipster.tools.path, "mothur-silva-reference", "silva.bacteria"))
	template.path <- c(file.path(data.path, "silva.gold.ng.fasta"))
	chimera.options <- paste(chimera.options, ref,"=", template.path, sep="")
}
if (reference=="full"){
	# whole reference
	data.path <- c(file.path(chipster.tools.path,"mothur-silva-reference", "mothur-silva-reference-whole"))
	template.path <- c(file.path(data.path, "silva.gold.align")) 
	chimera.options <- paste(chimera.options, ref,"=", template.path, sep="")
}
if (reference=="none"){
	if (file.exists("a.count_table")){
		chimera.options <- paste(chimera.options, " count=a.count_table", sep="")
	}
	else {
		stop('CHIPSTER-NOTE: In order to use the denovo option, note that you have to give a count table as input!')
	}
}
	

# if (file.exists("a.names")){
#	uchime.options <- paste(uchime.options, " name=a.names", sep=",")
#}
# if (file.exists("a.groups")){
# 	uchime.options <- paste(uchime.options, " group=a.groups", sep=",")
# }

if (dereplicate=="false"){
	chimera.options <- paste(chimera.options, ", dereplicate=F", sep="")
}
if (dereplicate=="true"){
	chimera.options <- paste(chimera.options, ", dereplicate=T", sep="")
}

chimera.options <- paste(chimera.options, ")", sep="")

#stop(paste('CHIPSTER-NOTE: ', chimera.options))

# Write batch file
write(chimera.options, "batch.mth", append=F)

# command
command <- paste(binary, "batch.mth", "> log.txt 2>&1")
# system("ls -l > log.txt")
# system("cat *.logfile >> log.txt")

# run
system(command)

#Output File Names: 
#a.ref.uchime.chimeras
#a.ref.uchime.accnos
#a.denovo.uchime.chimeras
#a.denovo.uchime.accnos

# rename outputs in a unified way for the subsequent steps
if (method=="vsearch"){
	system("mv a.denovo.vsearch.accnos denovoaccnosfile")
	system("mv a.ref.vsearch.accnos refaccnosfile")
}
if (method=="uchime"){
	system("mv a.denovo.uchime.accnos denovoaccnosfile")
	system("mv a.ref.uchime.accnos refaccnosfile")
}

# batch file 2 for Remove.seqs to remove chimeric sequences from the fasta file and the count file
# note that one should add parameter dups=F so that dereplicate=T from the previous step would take effect (a sequence that was assigned as chimeric in one sample would be removed only from that sample). This requires names file, count file doesn't seem to be supported yet.
if (reference=="none"){
	write("remove.seqs(accnos=denovoaccnosfile, fasta=a.fasta, count=a.count_table)", "remove.mth", append=F)
} else {
	if(file.exists("a.count_table")){
		write("remove.seqs(accnos=refaccnosfile, fasta=a.fasta, count=a.count_table)", "remove.mth", append=F)
	}else {
		write("remove.seqs(accnos=refaccnosfile, fasta=a.fasta)", "remove.mth", append=F)
	}
}

# command
command2 <- paste(binary, "remove.mth", "> log2.txt 2>&1")

# run
system(command2)

# Post process output, make Summary.seqs command
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