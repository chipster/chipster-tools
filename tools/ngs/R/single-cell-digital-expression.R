# TOOL single-cell-digital-expression.R: "Create digital gene expression matrix" (Corrects bead synthesis errors and extracts gene expression values from a BAM file where reads have been tagged with gene names. The resulting Digital Gene Expression matrix, DGE, can be used for further analysis with the Seural tools.)
# INPUT input.bam: "prepared BAM" TYPE GENERIC
# OUTPUT OPTIONAL digital_expression.txt.gz
# OUTPUT OPTIONAL digital_expression.tsv
# OUTPUT OPTIONAL digital_expression_summary.txt
# OUTPUT OPTIONAL cleaned.bam
# OUTPUT OPTIONAL synthesis_stats.txt
# OUTPUT OPTIONAL synthesis_stats_summary.txt
# PARAMETER OPTIONAL num.barcodes: "Number of barcodes to correct for bead synthesis error" TYPE INTEGER DEFAULT 2000 (Roughly 2x the expected number of cells. The number of barcodes on which to perform the correction. It is advisable to use roughly 2 times the anticipated number cells, as it was empirically found out that this allows to recover nearly every defective cell barcode that corresponds to a STAMP, rather than an empty bead cell barcode.)
# PARAMETER OPTIONAL primer.sequence: "Sequence" TYPE STRING DEFAULT AAGCAGTGGTATCAACGCAGAGTGAATGGG (Sequence to trim off. As a default, SMART adapter sequence.)
# PARAMETER OPTIONAL filtering.type: "How to filter the DGE matrix" TYPE [MIN_NUM_GENES_PER_CELL:"Min number of genes per cell" , NUM_CORE_BARCODES:"Number of core barcodes"] DEFAULT MIN_NUM_GENES_PER_CELL (How to filter the DGE matrix, based on minimum number of reads per cell, or by choosing the top N cells with most reads. Set the number in the Filtering parameter field below.)
# PARAMETER OPTIONAL filter.param: "Filtering threshold" TYPE INTEGER DEFAULT 0 (The corresponding parameter for filtering the DGE matrix.)


# OUTPUT OPTIONAL log.txt
# PARAMETER OPTIONAL num.core.barcodes: "Number of core barcodes" TYPE INTEGER DEFAULT 100 (How many reads per cell barcode required)
# PARAMETER OPTIONAL num.genes: "Number of genes per cell" TYPE INTEGER DEFAULT 0 (How many genes per cell required)
# PARAMETER OPTIONAL num.transcripts: "Number of transcripts per cell" TYPE INTEGER DEFAULT 0 (How many transcripts per cell required)
#  MIN_NUM_TRANSCRIPTS_PER_CELL:"Min number of transcripts per cell" 

# ML 12.10.2016 created
# ML 09.05.2017 combined detecting bead synthesis error here
# ML 04.07.2017 added num.transcripts parameter

path.dropseq <- c(file.path(chipster.tools.path, "drop-seq_tools"))

if (filter.param < 1) { 
	stop(paste('CHIPSTER-NOTE: ', "Please set a reasonable number to the filtering parameter."))	
}

# STEP 1: Detect bead synthesis errors:
# command 
command <- paste(path.dropseq, "/DetectBeadSynthesisErrors I=input.bam O=cleaned.bam OUTPUT_STATS=synthesis_stats.txt SUMMARY=synthesis_stats_summary.txt NUM_BARCODES=",num.barcodes," PRIMER_SEQUENCE=", primer.sequence, " 2>>log.txt", sep="")
# run the tool
system(command)


# STEP 2: Digital expression matrix:
# command start
command.start <- paste(path.dropseq, "/DigitalExpression I=cleaned.bam O=digital_expression.txt.gz SUMMARY=digital_expression_summary.txt", sep="")
# parameters
command.parameters <- paste(filtering.type, "=", filter.param, sep="")
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

