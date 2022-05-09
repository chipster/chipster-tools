# TOOL combine-tar.R: "Combine tar files" (Combines tar files into a single tar file.)
# INPUT tar{...}.tar: "Tar files" TYPE GENERIC
# OUTPUT combined.tar

source(file.path(chipster.common.path, "zip-utils.R"))

# check out if the file is compressed and if so unzip it
input.names <- read.table("chipster-inputs.tsv", header=F, sep="\t")
for (i in 1:nrow(input.names)) {
	unzipIfGZipFile(input.names[i,1])	
}

# Concatenat tar files. Files need to be concatenated pairwise or the output is corrupted.
for (i in 1:(nrow(input.names)-1)) {
  system(paste("tar --concatenate --file=", input.names[1,1], " ", input.names[i+1,1], sep="",collapse=""))
}

# Rename output
system(paste("mv",input.names[1,1],"combined.tar"))