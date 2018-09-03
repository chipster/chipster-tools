# TOOL mothur-chimera-uchime.R: "Remove chimeric sequences" (Removes chimeric sequences from a fasta-formatted alignment using the VSEARCH or the UCHIME method. You can use the 16S rRNA Silva gold bacterial set as a reference, or you can detect chimeras de novo using the more abundant sequences in your samples as a reference. This tool is based on the Mothur tools Chimera.vsearch, Chimera.uchime and Remove.seqs. Please note that it can take some time to run this tool!)
# INPUT a.fasta: "FASTA file" TYPE FASTA
# INPUT OPTIONAL a.count_table: "Count table" TYPE GENERIC
# OUTPUT OPTIONAL chimeras.removed.fasta
# OUTPUT OPTIONAL chimeras.removed.summary.tsv
# OUTPUT OPTIONAL chimeras.removed.count_table
# PARAMETER OPTIONAL reference: "Reference" TYPE [bacterial: "16S rRNA Silva gold bacteria", none: "none, de novo"] DEFAULT bacterial (You can use the 16S rRNA Silva gold bacterial set as a reference, or you can detect chimeras de novo using the more abundant sequences in your samples as a reference. Note that if you choose none, you have to give a count table as input.)
# PARAMETER OPTIONAL method: "Method" TYPE [vsearch, uchime] DEFAULT vsearch (Chimera detection method to use. Note that VSEARCH is much faster than the UCHIME method.)
# PARAMETER OPTIONAL dereplicate: "Dereplicate" TYPE [false, true] DEFAULT false (De novo chimera detection uses the more abundant sequences from the same sample to check the query sequence. When a sequence is flagged as chimeric in one sample, it can be removed from only that sample by setting dereplicate = true, or from all samples by setting dereplicate = false.)

# OUTPUT OPTIONAL log.txt
# OUTPUT OPTIONAL log2.txt

# EK 18.06.2013
# 30.3.2016 ML added input options and parameters
# ML 21.12.2016 update (new version, new Silva version)
# ML 14.3.2017 reference option (bacterial vs whole) + count-file output
# AMS 30.5.2017 added the possibility to use more processors
# EK 1.6.2017 added the vsearch method. Removed the input option for names and groups file as we use the more compact count file now for duplicates.
# EK 9.6.2017 changed both methods to use the fasta-formatted reference silva.gold.ng.fasta, because vsearch doesn't worked with the aligned format (silva.gold.align). Removed the full reference as it wasn't really full.
# EK 3.9.2018 updated silva.gold path

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

# start building the command for the Mothur chimera.uchime and chimera.vsearch tools based on the method chosen.
if (method=="vsearch"){
	chimera.options <- paste("chimera.vsearch(fasta=a.fasta, processors=", chipster.threads.max, ",", sep="")
}
if (method=="uchime"){
	chimera.options <- paste("chimera.uchime(fasta=a.fasta, processors=", chipster.threads.max, ",", sep="")
}

if (reference=="bacterial"){
	# bacterial reference in fasta format
	data.path <- c(file.path(chipster.tools.path, "mothur-silva-reference", "silva-gold"))
	reference.path <- c(file.path(data.path, "silva.gold.ng.fasta"))
	chimera.options <- paste(chimera.options, "reference=", reference.path, sep="")
}
# if (reference=="full"){
	# whole reference
#	data.path <- c(file.path(chipster.tools.path,"mothur-silva-reference", "mothur-silva-reference-whole"))
#	reference.path <- c(file.path(data.path, "silva.gold.align")) 
#	chimera.options <- paste(chimera.options, "reference=", reference.path, sep="")
#}

# count_table is used only when reference=none. Note that chimera.vsearch doesn't work if you give it both the reference and the count_table.
if (reference=="none"){
	if (file.exists("a.count_table")){
		chimera.options <- paste(chimera.options, " count=a.count_table", sep="")
	}
	else {
		stop('CHIPSTER-NOTE: Note that in order to use the denovo option, you have to give a count table as input!')
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
# According to the manual should add dups=F so that dereplicate=T from the previous step would take effect, but this seems to work correctly without (a sequence that was assigned as chimeric in one sample would be removed only from that sample). According to the manual this requires names file, but count file seems to work as well.
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