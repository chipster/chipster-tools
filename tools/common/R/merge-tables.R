# TOOL merge-tables.R: "Merge tables" (Merges tables using identifiers in the first column.)
# INPUT table{...}.tsv: "Tables to merge" TYPE GENERIC
# OUTPUT combined.tsv: "Merged table"
# PARAMETER OPTIONAL include.everything: "Include all the rows in the result file" TYPE [yes, no] DEFAULT no (Include also the non-matching lines in the result file.)

# JTT 22.10.2007
# EK 20.4.2015 clarified the script
# AMS 13.5.2016 modified to accept multiple inputs
# ML 18.4.2019 Modify to use tables-utils.R merge_tables function

source(file.path(chipster.common.path, "tables-utils.R"))

# Read the input names into a table
input.names <- read.table("chipster-inputs.tsv", header=F, sep="\t")
# Call merge_tables function in tables-utils.R :
merged <- merge_tables(input.names, include.everything)

# Writes out the combined table
write.table(merged, "combined.tsv", sep="\t", row.names=T, col.names=T, quote=F)

# EOF
