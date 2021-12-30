# TOOL mothur-seq-error.R: "Assess error rates" (Assess error rates if you have co-sequenced a mock community. This tool pulls out the sequences that were from the Mock sample specified with the group name parameter. Then it calculates the error rates and removes the Mock sequences from the fasta, count and taxonomy-assignment file. This tool is based on the Mothur tools get.groups, seq.error and remove.groups )
# INPUT a.fasta: "FASTA file" TYPE GENERIC
# INPUT a.count: "Count table" TYPE MOTHUR_COUNT
# INPUT reference.fasta: "Reference file" TYPE GENERIC
# INPUT taxonomy1.txt: "Taxonomy assignment" TYPE GENERIC
# OUTPUT mock-removed.count_table
# OUTPUT mock-removed.fasta.gz
# OUTPUT filtered-taxonomy-assignment.txt
# OUTPUT error-summary.tsv
# OUTPUT error-rate-summary.txt
# PARAMETER group: "The Mock group name" TYPE STRING DEFAULT "Mock"
# PARAMETER aligned: "Does the reference contain aligned sequences" TYPE [TRUE: yes, FALSE: no] DEFAULT FALSE (The aligned parameter allows you to specify whether your reference sequences are aligned. default=no)

# OUTPUT log.txt
# OUTPUT log2.txt
# OUTPUT log3.txt
# ES 03.08.2021
source(file.path(chipster.common.path,"tool-utils.R"))
source(file.path(chipster.common.path,"zip-utils.R"))

#check out if the file is compressed and if so unzip it
unzipIfGZipFile("a.fasta")

# binary
binary <- c(file.path(chipster.tools.path,"mothur","mothur"))
version <- system(paste(binary,"--version"),intern = TRUE)
documentVersion("Mothur",version)

# run get.groups
getgroups.options <- paste("get.groups(count=a.count, groups=",group,", fasta=a.fasta)",sep="")
# Write batch file
documentCommand(getgroups.options)
write(getgroups.options,"batch.mth",append = FALSE)
# command
command <- paste(binary,"batch.mth","> log.txt 2>&1")
runExternal(command)
system(command)

#Output File names: log.txt
#a.pick.count
#a.pick.fasta

# rename files
system("mv a.pick.count picked.count_table")
system("mv a.pick.fasta picked.fasta")

# run seq.error
seqerror.options <- paste("seq.error(fasta=picked.fasta, count=picked.count_table, reference=reference.fasta, aligned=",aligned,")",sep="")
documentCommand(seqerror.options)
write(seqerror.options,"batch2.mth",append = FALSE)
command2 <- paste(binary,"batch2.mth","> log2.txt 2>&1")
runExternal(command2)
system(command2)

#Output File Names: 
#a.pick.error.summary
#a.pick.error.seq
#a.pick.error.chimera
#a.pick.error.seq.forward
#a.pick.error.seq.reverse
#a.pick.error.count
#a.pick.error.matrix
#a.pick.error.ref

# run remove.groups
removegroups.options<- paste("remove.groups(count=a.count, fasta=a.fasta, taxonomy=taxonomy1.txt, groups=",group,")",sep="")
documentCommand(removegroups.options)
write(removegroups.options,"batch3.mth",append = FALSE)
command3 <- paste(binary,"batch3.mth","> log3.txt 2>&1")
runExternal(command3)
system(command3)

#Output File names: 
#a.pick.count
#a.pick.fasta
#taxonomy1.pick.txt

# Rename output files
system("mv a.pick.count mock-removed.count_table")
system("mv a.pick.fasta mock-removed.fasta")
system("mv taxonomy1.pick.txt filtered-taxonomy-assignment.txt")

# Post process output and combine output from get.groups and seq.error to error-rate-summary.txt
system("grep -A 1 Selected log.txt > summary1.txt")
system("grep -A 13 Multiply log2.txt > summary2.txt")
file.create("error-rate-summary.txt")
line <- paste("\nMothur get.groups command selected sequences based on the mock group name",group,":",sep=" ")
write(line,file="error-rate-summary.txt",append=TRUE)
summary_data1 <- readLines("summary1.txt")
for (row in summary_data1){
    write(row,file="error-rate-summary.txt",append=TRUE)
}
write("Error rate information:",file="error-rate-summary.txt",append=TRUE)
summary_data2 <- readLines("summary2.txt")
for (row in summary_data2){
    write(row,file="error-rate-summary.txt",append=TRUE)
}

#Run summary.seqs
summaryseq.options <- paste("summary.seqs(fasta=mock-removed.fasta, count=mock-removed.count_table")
summaryseq.options <- paste(summaryseq.options,", processors=",chipster.threads.max,")",sep = "")
documentCommand(summaryseq.options)
write(summaryseq.options,"summary.mth",append = FALSE)
command4 <- paste(binary,"summary.mth","> log_raw.txt")
system(command4)

# zip output fasta
system("gzip mock-removed.fasta")

# Post process output
system("grep -A 10 Start log_raw.txt > summary2.tsv")
# Remove one tab to get the column naming look nice:
system("sed 's/^		/	/' summary2.tsv > error-summary.tsv")
