# TOOL row_count.R: "Count rows" (Counts how many lines there are)
# INPUT input.txt: "Input file" TYPE GENERIC
# OUTPUT output.txt
# RUNTIME R-4.5.1
# TOOLS_BIN ""

# AO 15.5.2017

line <- system("wc -l input.txt", intern = TRUE)
write(line, file = "output.txt")
