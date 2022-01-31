# TOOL DEXseq.R: "Differential exon expression using DEXSeq" (Infers differential exon usage from RNA-seq data using the Bioconductor package DEXSeq. In order to prepare the input, run the tool \"Count aligned reads per exons for DEXSeq\" for your BAM files and combine the results to a count table using the tool \"Utilities - Define NGS experiment\". Please use the group column of the phenodata file to indicate your experimental groups. You need to have at least two biological replicates in each group.)
# INPUT countfile.tsv: "Count table" TYPE GENERIC 
# INPUT META phenodata.tsv: "Phenodata" TYPE GENERIC
# OUTPUT OPTIONAL dexseq-genes-with-significant-exons.tsv: dexseq-genes-with-significant-exons.tsv
# OUTPUT OPTIONAL dexseq-exons.pdf: dexseq-exons.pdf
# OUTPUT OPTIONAL dexseq-MAplot.pdf: dexseq-MAplot.pdf
# OUTPUT OPTIONAL dexseq-dispersion-plot.pdf: dexseq-dispersion-plot.pdf
# PARAMETER OPTIONAL organism: "Reference organism" TYPE ["FILES genomes/dexseq .DEXSeq.gtf"] DEFAULT "SYMLINK_TARGET genomes/dexseq/default .DEXSeq.gtf" (Which organism is your data from.)
# PARAMETER pvalue: "Threshold for adjusted p-value" TYPE DECIMAL FROM 0 TO 1 DEFAULT 0.05 (Threshold for BH adjusted p-values. If a gene has at least one exon below this p-value, all its exons will be included in the result list.)
# SLOTS 4

# 18.07.2013 JTT, Created
# 25.04.2014 MK, Modified for R-3.0
# 04.07.2014 AMS, New genome/gtf/index locations & names
# 08.07.2015 EK, Removed dexseq-all-genes.tsv from the results.
# 15.07.2016 ML, switched tools to newer versions

# Loads the library 
library(DEXSeq)

# Reads the phenodata
phenodata <- read.table("phenodata.tsv", header=T, sep="\t")
phenodata$condition<-phenodata$group
pd<-phenodata[,"group", drop=FALSE]
rownames(pd)<-phenodata$sample
colnames(pd)<-"condition"
pd$condition<-factor(pd$condition)

if(any(as.vector(table(phenodata$group))<2)) {
	stop("You need to have replicates for all groups you have specified.")
}

# Path to the gff file
gtf <- file.path(chipster.tools.path, "genomes", "dexseq", paste(organism, ".DEXSeq.gtf" ,sep="" ,collapse=""))

# Reads the data
d<-read.table("countfile.tsv", header=TRUE, sep="\t")
d2<-d[,grep("chip", colnames(d))]
cn<-substr(colnames(d2), 6, nchar(colnames(d2)))
for(i in 1:ncol(d2)) {
	v<-d2[,i, drop=F]
	rownames(v)<-rownames(d2)
	write.table(v, paste(cn[i], ".jtt", sep=""), col.names=FALSE, row.names=TRUE, sep="\t", quote=FALSE)
}


sampleTable = data.frame(
		row.names = phenodata$sample,
		condition = as.factor(phenodata$group))

ecs = DEXSeqDataSetFromHTSeq(countfiles = dir(pattern="jtt"), sampleData=sampleTable, design = ~ sample + exon + condition:exon, flattenedfile = gtf)


# Normalization
ecs<-estimateSizeFactors(ecs)

# Estimate dispersion

ecs<-estimateDispersions(ecs)

# Testing for differential exon usage
ecs <- testForDEU(ecs)
ecs = estimateExonFoldChanges(ecs, fitExpToVar="condition")
res <- DEXSeqResults(ecs)

siggenes<-as.character(unique(res$groupID[res$padj<pvalue]))
res2<-res[as.character(res$groupID) %in% siggenes,]

write.table(res2, "dexseq-genes-with-significant-exons.tsv", col.names=TRUE, row.names=TRUE, sep="\t", quote=FALSE)

# Visualization
# Note: as many pages of PDF as there are DEgenes.
if(nrow(res2)>0) {
	genes<-unique(as.character(res[which(res$padj<=pvalue),]$groupID))
	pdf("dexseq-exons.pdf", width=297/25.4, height=210/25.4)
	for(i in 1:length(genes)) {
		plottry<-try(plotDEXSeq(res, genes[i], displayTranscripts = FALSE, cex.axis = 1.2, cex = 1.3, lwd = 2, legend = TRUE))
		if(class(plottry)=="try-error") {
			plot(x=1, y=1, xlab="", ylab="", axes=F, type="")
			title(main=genes[i])
			text(x=1, y=1, "No results to plot for this gene")
		}
	}
	dev.off()
}
if(nrow(res)>0) {
	pdf("dexseq-MAplot.pdf", width=297/25.4, height=210/25.4)
	plotMA(res)
	dev.off()
}

#plotDispEsts = function( cds, ymin, linecol="#ff000080",
#  xlab = "mean of normalized counts", ylab = "dispersion",
#  log = "xy", cex = 0.45, ... )
#{
#  px = rowMeans( counts( cds, normalized=TRUE ) )
#  sel = (px>0)
#  px = px[sel]
#
#  py = fData(cds)$dispBeforeSharing[sel]
#  if(missing(ymin))
#      ymin = 10^floor(log10(min(py[py>0], na.rm=TRUE))-0.1)
#
#  plot(px, pmax(py, ymin), xlab=xlab, ylab=ylab,
#    log=log, pch=ifelse(py<ymin, 6, 16), cex=cex, ... )
#  xg = 10^seq( -.5, 5, length.out=100 )
#  fun = function(x) { cds@dispFitCoefs[1] + cds@dispFitCoefs[2] / x }
#  lines( xg, fun(xg), col=linecol, lwd=4)
#}

pdf("dexseq-dispersion-plot.pdf", width=297/25.4, height=210/25.4)
plotDispEsts(ecs, cex=0.2)
title(main="Dispersion plot")
legend(x="topright", legend="fitted dispersion", col="red", cex=1, pch="-")
dev.off()

