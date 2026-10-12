# TOOL make_tar.R: "Make a tar package" (Makes a tar package with selected files. Note that file names must be unique.)
# INPUT file{...}.tsv: "Files to include" TYPE GENERIC
# OUTPUT OPTIONAL chipster.tar
# PARAMETER name: "File name for tar package" TYPE STRING DEFAULT "chipster" (File name for the tar package. Ending .tar will be added to the name.)
# RUNTIME R-4.5.1
# TOOLS_BIN ""

source(file.path(chipster.common.lib.path, "tool-utils.R"))

# Read input names
input.names <- read.table("chipster-inputs.tsv", header = F, sep = "\t", quote = "", comment.char = "#", colClasses = "character", na.strings = character(0))

# Check that the file names don't point outside the job folder
safe_file_name(input.names[, 2])

# Check for duplicate file names
if (anyDuplicated(input.names[2])) {
    message <- paste("You have selected files with duplicated file names. File names must be unique. Please rename the files.")
    stop(paste("CHIPSTER-NOTE: ", message))
}

# Renamefiles to display names
for (i in 1:nrow(input.names)) {
    system(paste("mv --backup=numbered --suffix=. --", shQuote(input.names[i, 1]), shQuote(input.names[i, 2])))
}

# Tar
system("tar --exclude=\'chipster-inputs.tsv\' -cf chipster.tar -- *")


# Handle output names
#
# Define output name
if (nchar(name) > 0) {
    filename <- name
} else {
    filename <- "chipster"
}

filename <- paste(filename, ".tar", sep = "")


# Make a matrix of output names
outputnames <- matrix(NA, nrow = 1, ncol = 2)
outputnames[1, ] <- c("chipster.tar", filename)

# Write output definitions file
write_output_definitions(outputnames)
