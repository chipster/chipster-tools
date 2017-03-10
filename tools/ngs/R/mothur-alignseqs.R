# TOOL mothur-alignseqs.R: "Align sequences with Mothur" (Given a fasta file of 16S rRNA sequences, aligns them to the Silva reference set. Kmer searching with 8mers is followed by Needleman-Wunsch pairwise alignment, which penalizes the same amount for opening and extending a gap. This tool is based on the Mothur package. You can give a customized reference fasta. If you do so, make sure the files are correctly assigned in the parameters section! Please note that it can take some time to run this tool.)
# INPUT reads.fasta: "FASTA file" TYPE FASTA
# INPUT OPTIONAL reference.fasta: "custom reference FASTA file" TYPE FASTA
# OUTPUT OPTIONAL aligned.fasta
# OUTPUT aligned-summary.tsv



# EK 05.06.2013
# ML 21.12.2016 update (new Silva version)
# ML 4.1.2016 new, whole Silva reference
# OUTPUT log.txt

# check out if the file is compressed and if so unzip it
source(file.path(chipster.common.path, "zip-utils.R"))
unzipIfGZipFile("reads.fasta")

# binary
binary <- c(file.path(chipster.tools.path, "mothur", "mothur"))
# new bacteria references:
#data.path <- c(file.path(chipster.tools.path, "mothur-silva-reference"))
#template.path <- c(file.path(data.path, "silva.bacteria/silva.bacteria.fasta"))
# new whole references:
data.path <- c(file.path(chipster.tools.path, "mothur-silva-reference-whole"))
template.path <- c(file.path(data.path, "silva.nr_v123.align")) # or silva.seed_v123.align ?? https://mothur.org/wiki/Alignment_database

# batch file
if (file.exists("reference.fasta")){
	write(paste("align.seqs(fasta=reads.fasta, reference=reference.fasta)", sep=""), "batch.mth", append=F)
}else {
write(paste("align.seqs(fasta=reads.fasta, template=", template.path, ")", sep=""), "batch.mth", append=F)
}

# command
command <- paste(binary, "batch.mth", "> log.txt 2>&1")


# run
system(command)

system("mv reads.align aligned.fasta")

# batch file 2
write("summary.seqs(fasta=aligned.fasta)", "summary.mth", append=F)

# command
command2 <- paste(binary, "summary.mth", "> log_raw.txt")

# run
system(command2)

# Post process output
system("grep -A 9 Start log_raw.txt > aligned-summary2.tsv")
# Remove one tab to get the column naming look nice:
system("sed 's/^		/	/' aligned-summary2.tsv > aligned-summary.tsv")
