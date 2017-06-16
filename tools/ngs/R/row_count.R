# TOOL row_count.R: "Row count" (Counts how many lines there are)
# INPUT input.txt: "Input file" TYPE GENERIC 
# OUTPUT output.txt

# AO 15.5.2017

line<-system("wc -l input.txt", intern = TRUE)
write(line, file = "output.txt")
