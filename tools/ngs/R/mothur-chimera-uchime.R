# TOOL mothur-chimera-uchime.R: "Remove chimeric sequences with Mothur" (Remove chimeric sequences from a fasta-formatted alignment using the uchime method and the 16S rRNA Silva gold reference set. This tool is based on the Mothur package.)
# INPUT a.fasta: "FASTA file" TYPE FASTA
# OUTPUT OPTIONAL chimeras-removed.fasta
# OUTPUT OPTIONAL chimeras-removed-summary.tsv

# EK 18.06.2013
# ML 21.12.2016 update (new version, new Silva version)

# check out if the file is compressed and if so unzip it
source(file.path(chipster.common.path, "zip-utils.R"))
unzipIfGZipFile("a.fasta")

# binary
binary <- c(file.path(chipster.tools.path, "mothur", "mothur"))
data.path <- c(file.path(chipster.tools.path, "mothur-silva-reference"))
template.path <- c(file.path(data.path, "silva.bacteria/silva.gold.ng.fasta"))

# batch file
write(paste("chimera.uchime(fasta=a.fasta, template=", template.path, ")", sep=""), "batch.mth", append=F)

# command
command <- paste(binary, "batch.mth", "> log_raw.txt 2>&1")

# run
system(command)

# batch file 2
# note: output name was changed: a.uhcime.accnos => a.ref.uchime.accnos
write("remove.seqs(accnos=a.ref.uchime.accnos, fasta=a.fasta)", "remove.mth", append=F)

# command
command2 <- paste(binary, "remove.mth", ">> log_raw.txt 2>&1")

# run
system(command2)

# Post process output
system("mv a.pick.fasta chimeras-removed.fasta")

# system("grep -A 2 Removed log_raw.txt > log.txt")

# batch file 3
write("summary.seqs(fasta=chimeras-removed.fasta)", "summary.mth", append=F)

# command 3
command3 <- paste(binary, "summary.mth", ">> log_raw.txt")

# run
system(command3)

# Post process output
system("grep -A 9 Start log_raw.txt > chimeras-removed-summary2.tsv")
# Remove one tab to get the column naming look nice:
system("sed 's/^		/	/' chimeras-removed-summary2.tsv > chimeras-removed-summary.tsv")