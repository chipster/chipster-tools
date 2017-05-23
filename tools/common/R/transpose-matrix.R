# TOOL transpose-matrix.R: "Transpose matrix" (Given a matrix X returns the transpose of X. Returns same file type)
# INPUT input: input TYPE GENERIC
# OUTPUT output: output
# PARAMETER  header_text: "Header" TYPE [yes: Yes, no: No] DEFAULT yes (Is there a header in the input file)

# PARAMETER OPTIONAL file_format: "File format" TYPE [tsv: tsv, csv: csv] DEFAULT tsv (In which format is the input file, currently only tsv is supported)

#AO 22.5.2017

## Handle output names
# Source read_input_definitions and strip_name functions
source(file.path(chipster.common.path, "tool-utils.R"))
# read input names and strip file extension
input.names <- read_input_definitions()
input.name1 <- input.names$input

# Make a matrix of output names
# This overrides the default ones
output.names <- matrix(NA, nrow=1, ncol=2)
output.names[1,] <- c("output", input.name1)

# Write output definitions file
write_output_definitions(output.names)

##Check the header parameter
if( header_text == "yes") {
	header.bool <- TRUE
} else {
	header.bool <- FALSE
}

## Load, transpose, and write back
	
#if(file_format == "tsv") {
	matrix <- read.table(file = "input", sep = '\t', header = header.bool)
	matrix.transposed <- t(matrix)
	write.table(matrix.transposed, file = "output", sep = '\t', col.names = header.bool, row.names = header.bool)	
#} else if(file_format == "csv") {
	#matrix <- read.csv(input, HEADER = header.bool)
	#matrix.transposed <- t(matrix)
	#write.csv(matrix.transposed, file = "output", row.names=header.bool, col.names=header.bool)
#}
