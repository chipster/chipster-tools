# TOOL make_tar.R: "Make a tar package" (Makes a tar package with selected files. Note that file names must be unique.)
# INPUT file{...}.tsv: "Files to include" TYPE GENERIC
# OUTPUT OPTIONAL chipster.tar
# PARAMETER name: "File name for tar package" TYPE STRING DEFAULT "chipster" (File name for the tar package. Ending .tar will be added to the name.)
# IMAGE comp-20.04-r-deps
# RUNTIME R-4.1.1

# Read input names
input.names <- read.table("chipster-inputs.tsv", header=F, sep="\t")

# Check for duplicate file names
if (anyDuplicated(input.names[2])){
	message <- paste("You have selected files with duplicated file names. File names must be unique. Please rename the files.")
	stop(paste('CHIPSTER-NOTE: ', message ))
}


# Renamefiles to display names
for (i in 1:nrow(input.names)) {
	system(paste("mv --backup=numbered --suffix=.", input.names[i,1], input.names[i,2]))
}

# Tar
system("tar --exclude=\'chipster-inputs.tsv\' -cf chipster.tar *")


# Handle output names
#
source(file.path(chipster.common.path, "tool-utils.R"))

# Define output name
if (nchar(name) > 0){
	filename <- name
}else{
	filename <- "chipster"
}

filename <- paste(filename, ".tar", sep = "")


# Make a matrix of output names
outputnames <- matrix(NA, nrow=1, ncol=2)
outputnames[1,] <- c("chipster.tar", filename)

# Write output definitions file
write_output_definitions(outputnames)



