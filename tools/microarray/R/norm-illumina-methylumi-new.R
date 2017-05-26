# TOOL norm-illumina-methylumi-new.R: "Illumina - methylumi pipeline NEW" (Illumina methylation assay normalization using lumi methodology. As input, give simple tab delimited text file with the headers from Beadstudio, and a similarly constructed control file, OR a directly imported FinalReport file.  )
# INPUT samples.txt: samples.txt TYPE GENERIC 
# INPUT OPTIONAL controls.txt: controls.txt TYPE GENERIC 
# OUTPUT normalized.tsv: normalized.tsv 
# OUTPUT unmethylated.tsv: unmethylated.tsv 
# OUTPUT methylated.tsv: methylated.tsv 
# OUTPUT META phenodata.tsv: phenodata.tsv 
# OUTPUT OPTIONAL QC-plot.pdf: QC-plot.pdf
# PARAMETER chiptype: Chiptype TYPE [FDb.InfiniumMethylation.hg19:FDb.InfiniumMethylation.hg19] DEFAULT FDb.InfiniumMethylation.hg19 (Select the annotation package)
# PARAMETER OPTIONAL image.width: "Image width" TYPE INTEGER FROM 200 TO 3200 DEFAULT 600 (Width of the QC image)
# PARAMETER OPTIONAL image.height: "Image height" TYPE INTEGER FROM 200 TO 3200 DEFAULT 600 (Height of the QC image)
# PARAMETER OPTIONAL Pval: "Min. average p-val for filtering" TYPE DECIMAL FROM 0 TO 1 DEFAULT 0.05 (Remove poor quality samples, i.e. remove sample if the average detection p-val is above this threshold.)

# Old parameters:
# PARAMETER OPTIONAL normalization: Normalization TYPE [none: none, quantile: quantile, ssn: ssn] DEFAULT quantile ()
# PARAMETER OPTIONAL color.balance.adjustment: "Color balance adjustment" TYPE [none: none, quantile: quantile, ssn: ssn] DEFAULT quantile (Adjustment of color balance)
# PARAMETER OPTIONAL background.correction: "Background correction" TYPE [none: none, bgAdjust2C: bgAdjust2C, forcePositive: forcePositive] DEFAULT none (Should background adjustment be applied)
# PARAMETER OPTIONAL QCplots: QCplots TYPE [yes: yes, no: no] DEFAULT yes (Do you want quality control plots)
# PARAMETER target: "Illumina ID" TYPE STRING DEFAULT TargetID (Name of the TargetID column)
# PARAMETER signala: "Signal_A/Grn pattern" TYPE STRING DEFAULT Signal_A (Pattern identifying and common to all unmethylated data columns)
# PARAMETER signalb: "Signal_B/Red pattern" TYPE STRING DEFAULT Signal_B (Pattern identifying and common to all methylated data columns)


# ML: 17.5.2017: Update norm-illumina-methyllumi.R according to this: https://www.bioconductor.org/packages/devel/bioc/vignettes/methylumi/inst/doc/methylumi.pdf

# Loading libraries
library(lumi)
library(methylumi)
library(annotate)

if(chiptype=="FDb.InfiniumMethylation.hg19") {
	chiptype <- "FDb.InfiniumMethylation.hg19"
}

# Read data in
samples_file <- c("samples.txt")
if(file.exists("controls.txt")) {
	controls_file <- c("controls.txt")
	dat <- methylumiR(samples_file, qcfile=controls_file)
}else {
	dat <- methylumiR(samples_file)
}

# normalize:
avgPval <- colMeans(pvals(dat))
toKeep <- (avgPval<Pval)
# pData(mldat)$Gender[9] <- "F"  testidataspesifi homma
mldat.norm <- normalizeMethyLumiSet(dat[,toKeep])

dat4 <- mldat.norm

# QC plots
pdf(file="QC-plot.pdf", width=image.width/72, height=image.height/72)
par(mfrow=c(2,2))
plotColorBias1D(dat, main="Unpreprocessed")
#plotColorBias1D(dat4, main="Preprocessed")
plotColorBias1D(mldat.norm, main="Preprocessed")
boxplotColorBias(dat, main="Unpreprocessed")
# boxplotColorBias(dat4, main="Preprocessed")
boxplotColorBias(mldat.norm, main="Preprocessed")
par(las=2)
barplot(avgPval,ylab="Average P-Value")
dev.off()


# Convert sample names to Chipster style
# miksi tehty näin monta kertaa??
dat5<-exprs(dat4)
sample.names<-colnames(dat5)
sample.names<-paste("chip.", sample.names, sep="")
names(dat5)<-sample.names
colnames(dat5)<-sample.names

dat6<-methylated(dat4)
sample.names<-colnames(dat6)
sample.names<-paste("chip.", sample.names, sep="")
names(dat6)<-sample.names
colnames(dat6)<-sample.names

dat7<-unmethylated(dat4)
sample.names<-colnames(dat7)
sample.names<-paste("chip.", sample.names, sep="")
names(dat7)<-sample.names
colnames(dat7)<-sample.names

# annotations:
symbols <- mldat.norm@featureData@data$SYMBOL

# Write out a phenodata
group<-c(rep("", ncol(dat5)))
write.table(data.frame(sample=sample.names, chiptype=chiptype, group=group), file="phenodata.tsv", sep="\t", row.names=F, col.names=T, quote=F)

# Write out expression data
if(isEmpty(symbols)) {
	write.table(data.frame(dat5), file="normalized.tsv", col.names=T, quote=F, sep="\t", row.names=T)
	write.table(data.frame(dat6), file="methylated.tsv", col.names=T, quote=F, sep="\t", row.names=T)
	write.table(data.frame(dat7), file="unmethylated.tsv", col.names=T, quote=F, sep="\t", row.names=T)
} else {
	write.table(data.frame(symbol=symbols, dat5), file="normalized.tsv", col.names=T, quote=F, sep="\t", row.names=T)
	write.table(data.frame(symbol=symbols, dat6), file="methylated.tsv", col.names=T, quote=F, sep="\t", row.names=T)
	write.table(data.frame(symbol=symbols, dat7), file="unmethylated.tsv", col.names=T, quote=F, sep="\t", row.names=T)
	
}
