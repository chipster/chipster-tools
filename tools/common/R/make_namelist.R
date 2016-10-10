# TOOL make_namelist.R: "Make a list of file names" (Makes a list of file names of the selected files.)
# INPUT table{...}.file: "Input files" TYPE GENERIC
# OUTPUT files.txt 
# PARAMETER name: "File name for list" TYPE STRING DEFAULT "files.txt" (File name for the list.)
# PARAMETER sort: "Sort file" TYPE [yes, no] DEFAULT yes (Sort the list alphabetically.)

# Read input names
input.names <- read.table("chipster-inputs.tsv", header=F, sep="\t")

# Add original names to a vector
input.list <- vector(mode="character", length=0)
for (i in 1:nrow(input.names)) {
	input.list <- c(input.list, paste(input.names[i,2]))
}

# Sort
if (sort == "yes"){
	input.list <- sort(input.list)
}

# Write list
write.table(input.list, "files.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

# Handle output names
#
source(file.path(chipster.common.path, "tool-utils.R"))

# Define output name
if (nchar(name) > 0){
	filename <- name
}else{
	filename <- "files.txt"
}

# Make a matrix of output names
outputnames <- matrix(NA, nrow=1, ncol=2)
outputnames[1,] <- c("files.txt", filename)

# Write output definitions file
write_output_definitions(outputnames)



