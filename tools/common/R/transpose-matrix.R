# TOOL transpose-matrix.R: "Transpose matrix" (Given a matrix X returns the transpose of X. Returns same file type)
# INPUT input: input TYPE GENERIC
# OUTPUT output: output
# PARAMETER  header_text: "Does the first column have rownames" TYPE [yes: Yes, no: No] DEFAULT yes (Does the first column have rownames in the input file)
# PARAMETER OPTIONAL precision: "How many digits is used to express a number" DEFAULT 7 (How many digits is used to express a number)

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

## Check the header parameter
if( header_text == "yes") {
	header.bool <- TRUE
} else {
	header.bool <- FALSE
}

# Control number of digits, the default is 7
options(digits=precision)

## Load, transpose, and write back

matrix <- read.table(file = "input", sep = '\t', header = header.bool)
matrix.transposed <- t(matrix)
write.table(matrix.transposed, file = "output", sep = '\t', col.names = header.bool, row.names = header.bool, quote= FALSE)	
