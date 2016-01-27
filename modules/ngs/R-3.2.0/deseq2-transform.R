# TOOL deseq2-transform.R: "Transform read counts" (Transforms read counts for visualization and clustering purposes using the DESeq2 Bioconductor package. You can choose either vst or rlog transformation, or just DESeq2 normalisation. Note that the input file has to contain all the genes, not just differentially expressed ones. If you need to obtain transformed counts for differentially expressed genes, please click on \"More help\" to see the manual. You can create the input count table and phenodata file using the tool \"Utilities - Define NGS experiment\".)
# INPUT data.tsv: "Count table" TYPE GENERIC
# INPUT META phenodata.tsv: "Phenodata file" TYPE GENERIC
# OUTPUT OPTIONAL normalized-counts.tsv
# OUTPUT OPTIONAL vst-transformed-counts.tsv
# OUTPUT OPTIONAL rlog-transformed-counts.tsv
# PARAMETER method: "Transformation method" TYPE [vst: "variance stabilizing transformation (VST\)", rlog: "regularized log transformation (rlog\)", none: "no transformation, only DESeq2 normalization"] DEFAULT vst (Which method should be used to transform the data. VST runs faster than rlog. If the library size of the samples and therefore their size factors vary widely, the rlog transformation is a better option than VST.)

# EK 25.1.2016

# Load the library
library(DESeq2)

# Load the counts data and extract expression value columns
dat <- read.table("data.tsv", header=T, sep="\t", row.names=1)
dat2 <- dat[,grep("chip", names(dat))]

# Get the experimental group information from the phenodata. It is needed only for building the dds object and running the DESeq command, not for transformations.
phenodata <- read.table("phenodata.tsv", header=T, sep="\t")
condition <- as.character (phenodata[,pmatch("group",colnames(phenodata))])

# Create a DESeqDataSet object. The design is not needed, because the transformations are blind to the group information (blind=TRUE mode), meaning that dispersions are re-estimated using only an intercept.
dds <- DESeqDataSetFromMatrix(countData=dat2, colData=data.frame(condition), design = ~ condition)

# Run DESeq command to produce size factors and dispersions
dds <- DESeq(dds)

# no transformation
if (method == "none") {
	normalized.counts <- round(counts(dds,normalized=TRUE),digits=2)
	colnames(normalized.counts)<-colnames(dat2)
	write.table(normalized.counts,file="normalized-counts.tsv",sep="\t",row.names=T,quote=F)
	}

# VST transformation
if (method == "vst") {
	vst<-varianceStabilizingTransformation(dds)
	vstdf<-data.frame(round(assay(vst),digits=2))
	colnames(vstdf)<-colnames(dat2)
	write.table(vstdf,file="vst-transformed-counts.tsv",sep="\t",row.names=T,quote=F)
	}
	
# rlog transformation
if (method == "rlog") {
	rld<-rlog(dds)
	rldf<-data.frame(round(assay(rld),digits=4))
	colnames(rldf)<-colnames(dat2)
	write.table(rldf,file="rlog-transformed-counts.tsv",sep="\t",row.names=T,quote=F)
	}
	