# TOOL single-cell-digital-expression.R: "Digital expression" (Extracting Digital Gene Expression DGE data from an aligned library.)
# INPUT input.bam: "prepared BAM" TYPE GENERIC
# OUTPUT OPTIONAL digital_expression.txt.gz
# OUTPUT OPTIONAL digital_expression.tsv
# OUTPUT OPTIONAL digital_expression_summary.txt
# PARAMETER OPTIONAL num.barcodes: "Number of core barcodes" TYPE INTEGER DEFAULT 100 (How many barcodes)
# PARAMETER OPTIONAL num.genes: "Number of genes per cell" TYPE INTEGER DEFAULT 100 (How many genes per cell required)


# OUTPUT OPTIONAL log.txt


# ML 12.10.2016 created

path.dropseq <- c(file.path(chipster.tools.path, "drop-seq_tools"))

# command start
command.start <- paste(path.dropseq, "/DigitalExpression I=input.bam O=digital_expression.txt.gz SUMMARY=digital_expression_summary.txt", sep="")

# parameters
command.parameters <- paste("NUM_CORE_BARCODES=", num.barcodes, "MIN_NUM_GENES_PER_CELL=", num.genes)

# run the tool
command <- paste(command.start, command.parameters, " 2>> log.txt")
system(command)

# unzip the result file
system("gzip -d digital_expression.txt.gz 2>> log.txt")
system("mv digital_expression.txt digital_expression.tsv")

#digital_expression.tsv

# stop(paste('CHIPSTER-NOTE: ', command))

#EOF