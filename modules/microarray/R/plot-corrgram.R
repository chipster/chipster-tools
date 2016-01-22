# TOOL plot-corrgram.R: Correlogram (Plots the correlations between samples in a graph.)
# INPUT normalized.tsv: normalized.tsv TYPE GENE_EXPRS 
# INPUT META phenodata.tsv: phenodata.tsv TYPE GENERIC 
# OUTPUT OPTIONAL corrgram.png: corrgram.png 
# OUTPUT OPTIONAL corrgram.pdf: corrgram.pdf 
# PARAMETER OPTIONAL image.width: "Image width" TYPE INTEGER FROM 200 TO 3200 DEFAULT 600 (Width of the plotted network image)
# PARAMETER OPTIONAL image.height: "Image height" TYPE INTEGER FROM 200 TO 3200 DEFAULT 600 (Height of the plotted network image)
# PARAMETER OPTIONAL image.type: "Image type" TYPE [PNG: PNG, PDF: PDF] DEFAULT PNG (Format of the output image. Please note that if you choose the PDF format, it will take a lot of time to draw it, and thus it is advisable to export the PDF for closer examination.)


# Correlogram
# JTT 18.10.2007
# ML 22.1.2016, added option to print image either as png or pdf

# Parameter settings (default) for testing purposes
#image.width<-c(600)
#image.height<-c(600)

# Renaming variables
w<-image.width
h<-image.height

# Loads libraries
library(corrgram)

# Loads the normalized data
file<-c("normalized.tsv")
dat<-read.table(file, header=T, sep="\t", row.names=1)

# Separates expression values and flags
calls<-dat[,grep("flag", names(dat))]
dat2<-dat[,grep("chip", names(dat))]

# Loads phenodata
phenodata<-read.table("phenodata.tsv", header=T, sep="\t")

# Generating new labels
colnames(dat2)<-gsub(" ", "", phenodata$description)

# Plotting
if (image.type == "PNG") {
	bitmap(file="corrgram.png", width=w/72, height=h/72)
	corrgram(x=dat2, lower.panel=panel.shade, upper.panel=panel.pts, pch=".")
	dev.off()
}
if (image.type == "PDF") {
	pdf(file="corrgram.pdf", width=w/72, height=h/72)
	corrgram(x=dat2, lower.panel=panel.shade, upper.panel=panel.pts, pch=".")
	dev.off()
}



