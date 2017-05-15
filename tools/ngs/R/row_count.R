# TOOL axel_row_count_practice.R: "Row count" (Counts how many lines there is)
# INPUT input.txt: "Input file" TYPE GENERIC 
# OUTPUT output.txt

# AO 15.5.2017

line<-system("wc -l input.txt")
write(line, file = "output.txt")
