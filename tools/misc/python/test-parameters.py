# TOOL test-parameters.py: "Test parameters in Python" ()
# INPUT file{...}.tsv: "Files to include" TYPE GENERIC
# OUTPUT output.txt: ""
# PARAMETER OPTIONAL int TYPE INTEGER FROM -100 TO 100 DEFAULT 0
# PARAMETER OPTIONAL dec TYPE DECIMAL FROM -0.1 TO 1.0 DEFAULT 0
# PARAMETER OPTIONAL string: String TYPE STRING
# PARAMETER OPTIONAL enum: Enum TYPE [first: First, second: Second, third: Third] DEFAULT first
# PARAMETER OPTIONAL column: Column TYPE COLUMN_SEL
# PARAMETER OPTIONAL metacolumn: "Phenodata column" TYPE METACOLUMN_SEL
# PARAMETER OPTIONAL unchecked: "Unchecked string" TYPE UNCHECKED_STRING

# deprecated
# PARAMETER OPTIONAL file: File TYPE INPUT_SEL

with open('output.txt', 'w') as f:
    f.write('int	' + str(int) + '\n')
    f.write('dec	' + str(dec) + '\n')
    f.write('string	' + string + '\n')
    f.write('enum	' + str(enum) + '\n')
    f.write('column	' + str(column) + '\n')
    f.write('metacolumn	' + str(metacolumn) + '\n')
#	f.write('string	' + unchecked + '\n')
#	f.write('file	' + str(file) + '\n')
