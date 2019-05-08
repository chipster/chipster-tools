deseq_transform <- function(dat2, method){
	# Load the library
	library(DESeq2)
	
	# Get the experimental group information from the phenodata. It is needed only for building the dds object and running the DESeq command, not for transformations.
	phenodata <- read.table("phenodata.tsv", header=T, sep="\t")
	condition <- as.character (phenodata[,pmatch("group",colnames(phenodata))])
	
	# Create a DESeqDataSet object. The design is not needed, because the transformations are blind to the group information (blind=TRUE mode), meaning that dispersions are re-estimated using only an intercept.
	dds <- DESeqDataSetFromMatrix(countData=dat2, colData=data.frame(condition), design = ~ condition)
	
	# Run DESeq command to produce size factors and dispersions
	dds <- DESeq(dds)
	
	# no transformation
	if (method == "none") {
		transf.counts <- round(counts(dds,normalized=TRUE),digits=2)
		#colnames(normalized.counts)<-colnames(dat2)
		# write.table(normalized.counts,file="normalized-counts.tsv",sep="\t",row.names=T,quote=F)
		file.name <- "normalized-counts.tsv"
	}
	
	# VST transformation
	if (method == "vst") {
		vst<-varianceStabilizingTransformation(dds)
		transf.counts<-data.frame(round(assay(vst),digits=2))
		#colnames(vstdf)<-colnames(dat2)
		#write.table(vstdf,file="vst-transformed-counts.tsv",sep="\t",row.names=T,quote=F)
		file.name <- "vst-transformed-counts.tsv"
	}
	
	# rlog transformation
	if (method == "rlog") {
		rld<-rlog(dds)
		transf.counts<-data.frame(round(assay(rld),digits=4))
		#colnames(rldf)<-colnames(dat2)
		#write.table(rldf,file="rlog-transformed-counts.tsv",sep="\t",row.names=T,quote=F)
		file.name <- "rlog-transformed-counts.tsv"
	}
	
	colnames(transf.counts)<-colnames(dat2)
	return(transf.counts)
	write.table(transf.counts,file=file.name,sep="\t",row.names=T,quote=F)
	
}	