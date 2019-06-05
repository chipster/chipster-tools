# TOOL deseq2-transform.R: "Transform read counts" (Transforms read counts for visualization and clustering purposes using the DESeq2 Bioconductor package. You can choose either vst or rlog transformation, or just DESeq2 normalisation. Note that the input file has to contain all the genes, not just differentially expressed ones. If you need to obtain transformed counts for differentially expressed genes, please click on \"More help\" to see the manual. You can create the input count table and phenodata file using the tool \"Utilities - Define NGS experiment\".)
# INPUT data.tsv: "Count table" TYPE GENERIC
# INPUT META phenodata.tsv: "Phenodata file" TYPE GENERIC
# OUTPUT OPTIONAL transformed_counts.tsv
# PARAMETER method: "Transformation method" TYPE [vst: "variance stabilizing transformation (VST\)", rlog: "regularized log transformation (rlog\)", none: "no transformation, only DESeq2 normalization"] DEFAULT vst (Which method should be used to transform the data. VST runs faster than rlog. If the library size of the samples and therefore their size factors vary widely, the rlog transformation is a better option than VST.)

# EK 25.1.2016
# ML 8.5.2019 Move content to a separate function: also used by the heatmap tool


# Load the library
library(DESeq2)

# Load the counts data and extract expression value columns
dat <- read.table("data.tsv", header=T, sep="\t", row.names=1)
dat2 <- dat[,grep("chip", names(dat))]

# Call deseq_transform function:
source(file.path(chipster.common.path, "deseq-transform.R"))
transf.counts <- deseq_transform(dat2, method)

write.table(transf.counts, file="transformed_counts.tsv", sep="\t", row.names=T, quote=F)


# EOF
	