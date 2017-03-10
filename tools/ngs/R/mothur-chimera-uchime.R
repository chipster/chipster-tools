# TOOL mothur-chimera-uchime.R: "Remove chimeric sequences with Mothur" (Remove chimeric sequences from a fasta-formatted alignment using the uchime method and the 16S rRNA Silva gold reference set. This tool is based on the Mothur package. Please note that it can take some time to run this tool!)
# INPUT a.fasta: "FASTA file" TYPE FASTA
# INPUT OPTIONAL a.groups: "Groups file" TYPE MOTHUR_GROUPS
# INPUT OPTIONAL a.names: "Names file" TYPE MOTHUR_NAMES
# INPUT OPTIONAL a.count-table: "Count table" TYPE GENERIC
# OUTPUT OPTIONAL chimeras.removed.fasta
# OUTPUT OPTIONAL chimeras.removed.summary.tsv
# OUTPUT OPTIONAL chimeras.removed.count_table
# PARAMETER OPTIONAL dereplicate: "Dereplicate" TYPE [yes, no] DEFAULT yes (If sequence is flagged as chimeric, remove it from all samples)

# OUTPUT OPTIONAL log.txt
# OUTPUT OPTIONAL log2.txt

# EK 18.06.2013
# 30.3.2016 ML added input options and parameters
# ML 21.12.2016 update (new version, new Silva version)

# check out if the file is compressed and if so unzip it
source(file.path(chipster.common.path, "zip-utils.R"))
unzipIfGZipFile("a.fasta")

# binary
binary <- c(file.path(chipster.tools.path, "mothur", "mothur"))
# old references:
data.path <- c(file.path(chipster.tools.path, "mothur-data"))
template.path <- c(file.path(data.path, "silva.gold.align"))
# new bacterial references:
#data.path <- c(file.path(chipster.tools.path, "mothur-silva-reference"))
#template.path <- c(file.path(data.path, "silva.bacteria/silva.gold.ng.fasta"))
# new whole references:
#data.path <- c(file.path(chipster.tools.path, "mothur-data","mothur-silva-reference-whole"))
#template.path <- c(file.path(data.path, "silva.gold.align")) 


# batch file
uchime.options <- ""
uchime.options <- paste(uchime.options, "chimera.uchime(fasta=a.fasta, template=", template.path,sep="")
#uchime.options <- paste(uchime.options, "chimera.uchime(fasta=a.fasta", sep="")

if (file.exists("a.names")){
	uchime.options <- paste(uchime.options, " name=a.names", sep=",")
}
if (file.exists("a.groups")){
	uchime.options <- paste(uchime.options, " group=a.groups", sep=",")
}
if (file.exists("a.count-table")){
	uchime.options <- paste(uchime.options, " count=a.count_table", sep=",")
}
if (dereplicate=="yes"){
	uchime.options <- paste(uchime.options, ", dereplicate=F", sep="")
}
if (dereplicate=="no"){
	uchime.options <- paste(uchime.options, ", dereplicate=T", sep="")
}

uchime.options <- paste(uchime.options, ")", sep="")

# stop(paste('CHIPSTER-NOTE: ', uchime.options))

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

# batch file 2
write("remove.seqs(accnos=a.ref.uchime.accnos, fasta=a.fasta)", "remove.mth", append=F)

# command
command2 <- paste(binary, "remove.mth", "> log2.txt 2>&1")

# run
system(command2)

# Post process output
system("mv a.pick.fasta chimeras.removed.fasta")
if (file.exists("a.uchime.count_table")){
	system("mv a.uchime.count_table chimeras.removed.count_table")
}

# system("grep -A 2 Removed log_raw.txt > log.txt")

# batch file 3
write("summary.seqs(fasta=chimeras.removed.fasta)", "summary.mth", append=F)

# command 3
command3 <- paste(binary, "summary.mth", "> log_raw.txt")

# run
system(command3)

# Post process output
system("grep -A 10 Start log_raw.txt > chimeras.removed.summary2.tsv")
# Remove one tab to get the column naming look nice:
system("sed 's/^		/	/' chimeras.removed.summary2.tsv > chimeras.removed.summary.tsv")