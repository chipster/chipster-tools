# TOOL annotate-ensembl-ids.R: "Annotate Ensembl identifiers" (Annotates Ensembl IDs with gene symbols and descriptions, creates a new table containing these and the values in the original input file. The Ensembl IDs need to be either the rownames or in the first column of the input table.)
# INPUT genelist.tsv: genelist.tsv TYPE GENERIC
# OUTPUT annotated.tsv: annotated.tsv 
# PARAMETER species: Species TYPE [human: human, mouse: mouse, rat: rat] DEFAULT human (The species needs to be specified in order to annotate the Ensembl IDs.)


# ML, 05.02.2015
# ML, 12.11.2015, fixed
# ML, 14.08.2018, fixed + updated to new version
# JH (ML), 11.07.2019, fix entrezgene -> entrezgene_id
# ML, 06.11.2020, fix (runtime back to 3.2.3 from RUNTIME R-3.4.3)
# ML, 18.11.2020, add error message for ensembl connection issues

# ATM only for human, mouse, rat. 
# The ensembl IDs need to be either the row names or in the first column.

# Load libraries into memory
library(biomaRt)

# Load the data
file <- c("genelist.tsv")
#dat <- read.table(file, header=T, sep="\t", row.names=1)
dat <- read.table(file, header=T, sep="\t")

# Choose the ensembl IDs (whether in first column, or as rownames)
if(!is.na(pmatch("ENS", dat[2,1]))) {
	genes <- dat[,1]
}
if( !is.na(pmatch("ENS",  rownames(dat)[2] )) ||  !is.na(pmatch("ENS",  rownames(dat) )) ) {
	genes <- rownames(dat)
}
#if(!is.na(pmatch("ENS",  rownames(dat) ))) {
#	genes <- rownames(dat)
#}

if(is.na(pmatch("ENS", dat[2,1])) &&  is.na(pmatch("ENS",  rownames(dat) )) &&  is.na(pmatch("ENS",  rownames(dat)[2] )) ) {
	stop("CHIPSTER-NOTE: You can only annotate Ensembl IDs with this tool. The IDs need to be either rownames or in the first column of the table.")
}

# Fetch the gene symbols and descriptions from ENSEMBL using biomaRt
if (species=="human") {
	dataset <- "hsapiens_gene_ensembl"
	filt <- "hgnc_symbol" # not used?
}
if (species=="mouse") {
	dataset <- "mmusculus_gene_ensembl"
	filt <- "mgi_symbol"
}
if (species=="rat") {
	dataset <- "rnorvegicus_gene_ensembl"
	filt <- "rgd_symbol"
}

# If the connection can't be made for some reason, the error message is confusing; "Error in checkDataset(dataset = dataset, mart = mart): The given dataset:  xxxxx , is not valid.  Correct dataset names can be obtained with the listDatasets() function.""
# Check if this error comes up, and if so, ask user to try again later.
biomart_error_function <- function(e) {
	# This is the error one gets if there is connection problems:
  if(grepl("Error in checkDataset", e) ||  grepl("Error in useDataset", e)) {
	# Note: java expects to see string "Error: " in order to pass the CHIPSTER-NOTE. 
    print("Error: The original error message from biomaRt was: ")
    print(e)
	stop(paste("CHIPSTER-NOTE: ", "There seems to be a connection problem or lots of queries at the Ensembl mirrors. Please try again later."))
	# In case there is another, actual error, print that:
  }else { 
    print("There was problem with connecting to ensembl: ")
	print(e)
  }
}

# Retrieve the data from Ensembl (use the error function above to give nicer error to user):
# ensembl <- useMart("ensembl", dataset=dataset)
# ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset=dataset) # Let's leave this out and see if that helps: host="www.ensembl.org")
tryCatch(ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset=dataset), error = biomart_error_function ) #, finally = print("Hello"))

genes_ensembl_org <- getBM(attributes <- c("entrezgene_id", "ensembl_gene_id", "external_gene_name", "description"), filters = "ensembl_gene_id", values = genes, mart = ensembl, uniqueRows=T)

# Find the annotations for the identifiers in the input file:
pmatch_table		<- pmatch(genes, genes_ensembl_org[,2], duplicates.ok=T)
ensembl_table		<- as.data.frame(matrix(NA, nrow=length(genes), ncol=8))
ensembl_table[which(!is.na(pmatch_table)),] <- genes_ensembl_org[pmatch_table[(!is.na(pmatch_table))], ];
rownames(ensembl_table)	<- genes;
colnames(ensembl_table) <- colnames(genes_ensembl_org);

# Build the table:
# if identifiers in the first column:
if(!is.na(pmatch("ENS", dat[2,1]))) {
	results <- cbind(dat[,1], ensembl_table[,3:4], dat[,2:ncol(dat)]);
	colnames(results) <- c(colnames(dat)[1], "symbol", "description", colnames(dat)[2:ncol(dat)])
	# write result table to output
	write.table(results, file="annotated.tsv", col.names=T, quote=F, sep="\t", row.names=F)
}
# if identifiers = rownames:
if(!is.na(pmatch("ENS",  rownames(dat)[2] )) ||  !is.na(pmatch("ENS",  rownames(dat) )) ) {
	#results <- cbind(genes, ensembl_table[,3:4], dat);
	results <- cbind(ensembl_table[,3:4], dat);
	colnames(results) <- c("symbol", "description",  colnames(dat));
	rownames(results) <- rownames(dat)
	# write result table to output
	write.table(results, file="annotated.tsv", col.names=T, quote=F, sep="\t", row.names=T)
}




