# TOOL test-parameters.R: "Test parameters in R" ()  
# PARAMETER OPTIONAL int: Integer TYPE INTEGER FROM -100 TO 100 DEFAULT 0 
# PARAMETER OPTIONAL dec: Decimal TYPE DECIMAL FROM -0.1 TO 1.0 DEFAULT 0 
# PARAMETER OPTIONAL string: String TYPE STRING
# PARAMETER OPTIONAL unchecked: "Unchecked string" TYPE UNCHECKED_STRING 
# PARAMETER OPTIONAL enum: Enum TYPE [first: First, second: Second, third: Third] DEFAULT first 
# PARAMETER OPTIONAL column: Column TYPE COLUMN_SEL 
# PARAMETER OPTIONAL metacolumn: "Phenodata column" TYPE METACOLUMN_SEL 
# PARAMETER OPTIONAL file: File TYPE INPUT_SEL
