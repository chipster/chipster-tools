# TOOL deseq2-normalized-counts-for-gene.R: "Plot normalized counts for a gene" (Plot normalized counts for a gene \/genes using the DESeq2 Bioconductor package. You can type the gene name\(s\) as parameter, or give a list of names as a tsv file.)
# INPUT data.tsv: "Count table" TYPE GENERIC (Raw count data used for plotting.)
# INPUT OPTIONAL genelist.tsv: "List of genes for plotting" TYPE GENERIC (List of gene names for plotting as a tsv file. Note that giving this file as an input overwrites the list of genes given as a parameter!)
# INPUT META phenodata.tsv: "Phenodata file" TYPE GENERIC
# OUTPUT OPTIONAL normalized_counts.pdf
# PARAMETER OPTIONAL gene.names: "Gene names" TYPE STRING (Gene name\(s\). If you list multiple gene names, separate them with comma \(,\). Note that this list is NOT used if you give a list as a tsv file as input!)
# PARAMETER OPTIONAL show.names: "Show names in plot" TYPE [yes, no] DEFAULT yes (Show sample names in plot. In more complex cases this may make the plot too cluttered. To plot sample names, you need to determine them in the description column in phenodata file!)
# PARAMETER OPTIONAL how.many: "How many of the top genes from the input file are plotted" TYPE INTEGER FROM 1 DEFAULT 5 (If you give the genes to plot as an tsv file, this parameter sets the number of genes from the top of the table you wish to plot. Note that the pdf grows very large when you add more plots to it.)


# AMS 21.4.2015 
# ML 25.2.2021 Add option to give multiple gene names as input (as a list or tsv file)

# Loads the libraries
library(DESeq2)
library(ggplot2)

# Load the count table and extract expression value columns
dat <- read.table("data.tsv", header=T, sep="\t", row.names=1)
dat2 <- dat[,grep("chip", names(dat))]

# Get the experimental group information from the phenodata
phenodata <- read.table("phenodata.tsv", header=T, sep="\t")
condition <- as.character (phenodata[,pmatch("group",colnames(phenodata))])

# If a tsv file for gene names is given, use that. 

if (file.exists("genelist.tsv")) { 
  genelist <- read.table("genelist.tsv", header=T, sep="\t", row.names=1)
  # Choose first how.many genes or max number of rows of the input table for plotting
  if (nrow(genelist) < 5 ){ 
	  how.many <- nrow(genelist)
	}
  gene.names.clean <- rownames(genelist)[1:how.many]

}else{ 
# Otherwise, use the list of genes given as a parameter.
# Handle the input gene names
# Split from comma:
gene.names.split <- strsplit(gene.names, ",")[[1]]
# remove whitespace:
gene.names.clean <- gsub("[[:blank:]]", "", gene.names.split)
}

# Create a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData=dat2, colData=data.frame(condition), design = ~ condition)

dds <- DESeq(dds)
#res <- results(dds)
# desc <- phenodata[,7]

# Make plot as pdf
pdf(file="normalized_counts.pdf")
# One plot per gene:
for (i in 1:length(gene.names.clean) ) { 
	gene.name.to.print <- gene.names.clean[i]
  	d <- plotCounts(dds, gene=gene.name.to.print, intgroup="condition", returnData=TRUE)
  	if (show.names == "yes"){
		  # Check that there are descriptions in phenodata file that can be used in plotting:
			if (is.na(pmatch("description", colnames(phenodata)))) { 
  			stop("CHIPSTER-NOTE: To plot sample names, you need to determine them in the description column in phenodata file!")
			}else{
		  	desc <- phenodata[,"description"] 
			# Plotting:
   			# Note: as we are using ggplot within a for loop, we need to use "print"! 
    		print(ggplot(d, aes(x=condition, y=count, log="y")) +
    		#geom_point(color="blue", size=3, shape=5, position=position_jitter(w=0.1,h=0)) +
    		geom_point(color="blue", size=3, shape=5) +
    		geom_text(aes(label=desc),hjust=-0.5, vjust=0, color="black", size=4) +
    		ylab("normalized counts") +
    		xlab("group") +
    		ggtitle(gene.name.to.print)	)
			} 
  	}else{
		print(ggplot(d, aes(x=condition, y=count, log="y")) +
		#geom_point(color="blue", size=3, shape=5, position=position_jitter(w=0.1,h=0)) +
		geom_point(color="blue", size=3, shape=5) +
		ylab("normalized counts") +
		xlab("group") +
		ggtitle(gene.name.to.print) )
		}
} 
dev.off()
  



