# TOOL deseq2-pca-heatmap.R: "PCA and heatmap of samples with DESeq2" (Creates PCA and heatmap plots for samples using the DESeq2 Bioconductor package. Visualizing similarities and dissimilarities between samples can help you to perform experiment level quality control. This tool takes as input a table of raw counts. The count table has to be associated with a phenodata file describing the experimental groups. You can create the count table and phenodata file using the tool \"Utilities - Define NGS experiment\". You can distinguish the samples belonging to different sample groups based on different colors and/or shapes.)
# INPUT data.tsv: "Count table" TYPE GENERIC
# INPUT META phenodata.tsv: "Phenodata file" TYPE GENERIC
# OUTPUT OPTIONAL PCA_and_heatmap_deseq2.pdf
# PARAMETER column: "Phenodata column for coloring samples in PCA plot" TYPE METACOLUMN_SEL DEFAULT group (Phenodata column by which the samples will be colored in the PCA plot.)
# PARAMETER OPTIONAL second.column: "Phenodata column for the shape of samples in PCA plot" TYPE METACOLUMN_SEL DEFAULT EMPTY (Phenodata column by which the samples will be shaped in the PCA plot.)
# PARAMETER OPTIONAL show.names: "Show sample names in PCA plot" TYPE [yes, no] DEFAULT yes (You can add sample names to the PCA plot from the description column of the phenodata file. In more complex cases this may make the plot too cluttered.)

# EK 3.2.2015 
# AMS 22.4.2015 Added option for sample names in plot
# ML 4.1.2016 Added option to draw the PCA using two sample groups

# Loads the libraries
library(DESeq2)
library(gplots)
library(RColorBrewer)
library(ggplot2)

# Load the count table and extract expression value columns
dat <- read.table("data.tsv", header=T, sep="\t", row.names=1)
dat2 <- dat[,grep("chip", names(dat))]

# Get the experimental group information from the phenodata
phenodata <- read.table("phenodata.tsv", header=T, sep="\t")
condition <- as.character (phenodata[,pmatch(column,colnames(phenodata))])

if (second.column=="EMPTY") {
	dds <- DESeqDataSetFromMatrix(countData=dat2, colData=data.frame(condition), design = ~condition)
}
if (second.column!="EMPTY") {
	# condition2 <- as.character (phenodata[,pmatch(second.column,colnames(phenodata))])
	condition2 <- as.character (phenodata[,second.column])
	design <- data.frame(condition=as.factor(condition), condition2=as.factor(condition2))  
	rownames(design) <- colnames(dat2)
	#colnames(design) <- c(column, second.column)
	dds <- DESeqDataSetFromMatrix(countData=dat2, colData=design, design = ~ condition2 + condition)
	
}


# Create a DESeqDataSet object
#dds <- DESeqDataSetFromMatrix(countData=dat2, colData=data.frame(condition), design = ~ condition)

# Calculate size factors and estimate dispersions
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)

# Perform transformation
vst<-varianceStabilizingTransformation(dds)
vstmat<-assay(vst)

# Make PCA plot as pdf

# one group
if (second.column=="EMPTY") {
	data <- plotPCA(vst, intgroup=c("condition"), returnData=TRUE)
}
# two groups
if (second.column!="EMPTY") {
	data <- plotPCA(vst, intgroup=c("condition", "condition2"), returnData=TRUE)
	percentVar <- round(100 * attr(data, "percentVar"))
	desc <- phenodata[,"description"]
	ggplot(data, aes(PC1, PC2, color=condition, shape=condition2, title="PCA plot for TREATMENT and READTYPE")) +
			geom_point(size=6) +
			geom_text(aes(label=desc),hjust=0, vjust=1.7, color="black", size=4) +
			xlab(paste0("PC1: ",percentVar[1],"% variance")) +
			ylab(paste0("PC2: ",percentVar[2],"% variance")
			)
}

percentVar <- round(100 * attr(data, "percentVar"))
desc <- phenodata[,"description"]

pdf(file="01-pca-deseq2.pdf")
if (second.column=="EMPTY") {
	if (show.names == "yes"){
		ggplot(data, aes(PC1, PC2, color=condition)) +
				geom_point(size=6) +
				geom_text(aes(label=desc),hjust=0, vjust=1.7, color="black", size=4) +
				xlab(paste0("PC1: ",percentVar[1],"% variance")) +
				ylab(paste0("PC2: ",percentVar[2],"% variance"))
	}else{
		ggplot(data, aes(PC1, PC2, color=condition)) +
				geom_point(size=6) +
				xlab(paste0("PC1: ",percentVar[1],"% variance")) +
				ylab(paste0("PC2: ",percentVar[2],"% variance"))
	}
}
if (second.column!="EMPTY") {
	if (show.names == "yes"){
	ggplot(data, aes(PC1, PC2, color=condition, shape=condition2)) +
		geom_point(size=6) +
		geom_text(aes(label=desc),hjust=0, vjust=1.7, color="black", size=4) +
		xlab(paste0("PC1: ",percentVar[1],"% variance")) +
		ylab(paste0("PC2: ",percentVar[2],"% variance"))
}else{
	ggplot(data, aes(PC1, PC2, color=condition, shape=condition2)) +
			geom_point(size=6) +
			xlab(paste0("PC1: ",percentVar[1],"% variance")) +
			ylab(paste0("PC2: ",percentVar[2],"% variance"))
	}
}
dev.off()


# Make a distance matrix and name samples accroding to the phenodata description column
distvst<-dist(t(vstmat))
mdistvst<-as.matrix(distvst)
rownames(mdistvst)<-colnames(mdistvst)<-as.vector(phenodata$description)

# Calculate hierarchical clustering
hcvst<-hclust(distvst)

# Make the colors and plot the heatmap as pdf
hmcol<-colorRampPalette(brewer.pal(9,"GnBu"))(100)
pdf(file="02-heatmap-deseq2.pdf")
heatmap.2(mdistvst,Rowv=as.dendrogram(hcvst),Colv=as.dendrogram(hcvst),symm=TRUE,trace="none",col=rev(hmcol),margin=c(13,13))
dev.off()

# To be added: rlog transformation
# rld<-rlog(dds)
# rldmat<-assay(rld)
# distrld<-dist(t(rldmat))
# mdistrld<-as.matrix(distrld)
# hcrld<-hclust(distrld)
system("gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=PCA_and_heatmap_deseq2.pdf *.pdf")

