# TOOL test-parameters.R: "Test parameters in R" (Tool for testing all parameter types)
# INPUT normalized.tsv: "Input file" TYPE GENERIC (Input file description.)
# INPUT OPTIONAL normalized2.tsv: "Another optional file" TYPE GENERIC (Another input file description.)
# OUTPUT output.tsv: "Output file"
# PARAMETER OPTIONAL int: Integer TYPE INTEGER FROM -100 TO 100 DEFAULT 0 (integer parameter)
# PARAMETER OPTIONAL dec: Decimal TYPE DECIMAL FROM -0.1 TO 1.0 DEFAULT 0 (decimal parameter)
# PARAMETER OPTIONAL string: String TYPE STRING (string parameter)
# PARAMETER OPTIONAL unchecked: "Unchecked string" TYPE UNCHECKED_STRING (string parameter with special characters)
# PARAMETER OPTIONAL enum: Enum TYPE [first: First, second: Second, third: Third] DEFAULT first (multiple selection)
# PARAMETER OPTIONAL column: Column TYPE COLUMN_SEL (select column from input file)
# PARAMETER OPTIONAL metacolumn: "Phenodata column" TYPE METACOLUMN_SEL (select column from phenodata)
# PARAMETER OPTIONAL file: File TYPE INPUT_SEL (select one of the inputs)

print(paste("int: ", int))
print(paste("dec: ", dec))
print(paste("string: ", string))
print(paste("unchecked:", unchecked))
print(paste("enum:", enum))
print(paste("column:", column))
print(paste("metacolumn:", metacolumn))
print(paste("file:", file))

system("mv normalized.tsv output.tsv")
