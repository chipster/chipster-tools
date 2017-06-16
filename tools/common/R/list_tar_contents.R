# TOOL list_tar_contents.R: "List contents of a tar file" (List the contents of a tar package. The file can be gzip compressed.)
# INPUT file.tar: "Tar file" TYPE GENERIC
# OUTPUT toc.txt 

# List contents of tar to a file
system("tar tf file.tar > toc.txt")

# Change display name
source(file.path(chipster.common.path, "tool-utils.R"))
inputnames <- read_input_definitions()

basename <- strip_name(inputnames$file.tar)

# Make a matrix of output names
outputnames <- matrix(NA, nrow=1, ncol=2)
outputnames[1,] <- c("toc.txt", paste(basename, "_toc.txt", sep =""))

# Write output definitions file
write_output_definitions(outputnames)



