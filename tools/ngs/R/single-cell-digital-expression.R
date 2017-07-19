# TOOL single-cell-digital-expression.R: "Create digital gene expression matrix" (Extracting Digital Gene Expression DGE data from an aligned library.)
# INPUT input.bam: "prepared BAM" TYPE GENERIC
# OUTPUT OPTIONAL digital_expression.txt.gz
# OUTPUT OPTIONAL digital_expression.tsv
# OUTPUT OPTIONAL digital_expression_summary.txt
# OUTPUT OPTIONAL cleaned.bam
# OUTPUT OPTIONAL synthesis_stats.txt
# OUTPUT OPTIONAL synthesis_stats_summary.txt
# OUTPUT OPTIONAL log.txt
# PARAMETER OPTIONAL num.barcodes: "Number of barcodes" TYPE INTEGER DEFAULT 2000 (Roughly 2x the number of cells)
# PARAMETER OPTIONAL primer.sequence: "Sequence" TYPE STRING DEFAULT AAGCAGTGGTATCAACGCAGAGTGAATGGG (Sequence to trim off. As a default, SMART adapter sequence.)
# PARAMETER OPTIONAL num.core.barcodes: "Number of core barcodes" TYPE INTEGER DEFAULT 100 (How many barcodes)
# PARAMETER OPTIONAL num.genes: "Number of genes per cell" TYPE INTEGER DEFAULT 0 (How many genes per cell required)
# PARAMETER OPTIONAL num.transcripts: "Number of transcripts per cell" TYPE INTEGER DEFAULT 0 (How many transcripts per cell required)

# OUTPUT OPTIONAL log.txt


# ML 12.10.2016 created
# ML 09.05.2017 combined detecting bead synthesis error here
# ML 04.07.2017 added num.transcripts parameter

path.dropseq <- c(file.path(chipster.tools.path, "drop-seq_tools"))

# STEP 1: Detect bead synthesis errors:
# command 
command <- paste(path.dropseq, "/DetectBeadSynthesisErrors I=input.bam O=cleaned.bam OUTPUT_STATS=synthesis_stats.txt SUMMARY=synthesis_stats_summary.txt NUM_BARCODES=",num.barcodes," PRIMER_SEQUENCE=", primer.sequence, " 2>>log.txt", sep="")
# run the tool
system(command)


# STEP 2: Digital expression matrix:
# command start
command.start <- paste(path.dropseq, "/DigitalExpression I=cleaned.bam O=digital_expression.txt.gz SUMMARY=digital_expression_summary.txt", sep="")
# parameters
command.parameters <- paste("NUM_CORE_BARCODES=", num.core.barcodes, "MIN_NUM_GENES_PER_CELL=", num.genes, "MIN_NUM_TRANSCRIPTS_PER_CELL=", num.transcripts)
#command.parameters <- ""
#if (num.core.barcodes != "empty"){
#	command.parameters <- paste(command.parameters, "NUM_CORE_BARCODES=", num.core.barcodes)
#}else if (num.genes != "empty"){
#	command.parameters <- paste(command.parameters, "MIN_NUM_GENES_PER_CELL=", num.genes)
#}else if (num.transcripts != "empty"){
#	command.parameters <- paste(command.parameters, "MIN_NUM_TRANSCRIPTS_PER_CELL=", num.transcripts)
#}


# run the tool
command <- paste(command.start, command.parameters, " 2>> log.txt")
system(command)

# unzip the result file
system("gzip -d digital_expression.txt.gz 2>> log.txt")
system("mv digital_expression.txt digital_expression.tsv")

#digital_expression.tsv

# stop(paste('CHIPSTER-NOTE: ', command)


#EOF

