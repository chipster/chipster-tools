# TOOL gzip.R: "Compress a file with gzip" (Compress a file with gzip. Note that some file types like BAM are already compressed, and compressing them will not make the files any smaller.)
# INPUT file: "File to compress" TYPE GENERIC
# OUTPUT OPTIONAL file.gz 

source(file.path(chipster.common.path, "tool-utils.R"))

# We need to change the file name before compressing so it is preserved when uncompressing
input.names <- read.table("chipster-inputs.tsv", header=F, sep="\t")
system(paste("mv", input.names[1,1], input.names[1,2]))

# Gzip
system(paste("gzip", input.names[1,2]))

# Now we need to rename the gz file so that Chipster recognizes it as the output
filename <- paste(input.names[1,2], ".gz", sep = "")
system(paste("mv", filename, "file.gz"))

# And finally change Chipster display name to match original file name
outputnames <- matrix(NA, nrow=1, ncol=2)
outputnames[1,] <- c("file.gz", filename)

# Write output definitions file
write_output_definitions(outputnames)



