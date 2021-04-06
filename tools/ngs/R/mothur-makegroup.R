# TOOL mothur-makegroup.R: "Combine FASTA files and make a groups file" (Combines FASTA files of each sample to one file and creates a groups file. FASTA files can be gzipped. This tool is based on the Mothur tools make.group and merge.files.)
# INPUT reads{...}.fasta: "FASTA files" TYPE FASTA
# OUTPUT sequences.fasta.gz
# OUTPUT sequences.groups
# OUTPUT log1.txt

# EK 1.4.2021


source(file.path(chipster.common.path, "tool-utils.R"))
source(file.path(chipster.common.path, "zip-utils.R"))

# check out if the file is compressed and if so unzip it
unzipIfGZipFile("reads.fasta")

# Binary
binary <- c(file.path(chipster.tools.path,"mothur","mothur"))
version <- system(paste(binary,"--version"),intern = TRUE)
documentVersion("Mothur",version)

# Collect names of the input files
input.names <- read.table("chipster-inputs.tsv", header=F, sep="\t")

# Set names for first input
fasta_names <- input.names[1,1]
group_names <- strip_name(paste(input.names[1,2]))

# Loop through the rest of the inputs
for (i in 2:nrow(input.names)) {
  fasta_names <- paste(fasta_names, input.names[i,1], sep="-") 
  group <- strip_name(paste(input.names[i,2]))
  group_names <- paste(group_names, group, sep="-")
}

makegroup.options <- paste("make.group(fasta=",fasta_names, ", groups=",group_names, ")", sep = "")
documentCommand(makegroup.options)
write(makegroup.options,"makegroup.mth",append = FALSE)
command1 <- paste(binary,"makegroup.mth","> log1.txt")
system(command1)

# Rename groups file
filename <- list.files(path =".", pattern="*.groups")
system(paste("mv", filename[1], "sequences.groups"))

# Merge fasta files
mergefiles.options <- paste("merge.files(input=", fasta_names, ", output=sequences.fasta)",  sep = "")
documentCommand(mergefiles.options)
write(mergefiles.options,"mergefiles.mth",append = FALSE)
command2 <- paste(binary,"mergefiles.mth","> log2.txt")
system(command2)

# zip output fasta
system("gzip sequences.fasta")
#Sys.sleep(120)
