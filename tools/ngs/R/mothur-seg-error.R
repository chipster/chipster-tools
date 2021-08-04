# TOOL mothur-seq-error.R: "Assessing error rates"
# INPUT a.fasta: "FASTA file" TYPE GENERIC
# INPUT a.count: "Count table" TYPE GENERIC
# OUTPUT error.fasta.gz
# OUTPUT error.summary.tsv
# PARAMETER group: "The Mock group name" TYPE STRING DEFAULT Mock

# ES 03.08.2021

source(file.path(chipster.common.path,"tool-utils.R"))
source(file.path(chipster.common.path,"zip-utils.R"))

#check out if the file is compressed and if so unzip it
unzipIfGZipFile("a.fasta")
k
# binary
binary <- c(file.path(chipster.tools.path,"mothur","mothur"))
version <- system(paste(binary,"--version"),intern = TRUE)
documentVersion("Mothur",version)

getgroups.options <- paste("get.groups(count=a.count, fasta=a.fasta, groups=Mock)")
# Write batch file
documentCommand(getgroups.options)
write(getgroups.options,"batch.mth",append = FALSE)

# command
command <- paste(binary,"batch.mth","> log.txt 2>&1")
system(command)