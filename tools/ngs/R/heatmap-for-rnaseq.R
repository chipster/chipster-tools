# TOOL heatmap-for-rnaseq.R: "Heatmap for RNA-seq results" (Draw a heatmap and dendrogram of the genes of interest. As input, give the list of genes of interest and the original count table containing RAW COUNTS for all genes in all samples. The tool will first perform a transformation for the values.)
# INPUT de-list.tsv: "List of differentially expressed genes" TYPE GENERIC
# INPUT data.tsv: "Raw counts" TYPE GENERIC
# INPUT META phenodata.tsv: "Phenodata file" TYPE GENERIC 
# OUTPUT OPTIONAL log.txt
# OUTPUT OPTIONAL heatmap.pdf
# OUTPUT OPTIONAL vst-transformed-counts.tsv 
# PARAMETER method: "Transformation method" TYPE [vst: "variance stabilizing transformation (VST\)", rlog: "regularized log transformation (rlog\)", none: "no transformation, only DESeq2 normalization"] DEFAULT vst (Which method should be used to transform the data. VST runs faster than rlog. If the library size of the samples and therefore their size factors vary widely, the rlog transformation is a better option than VST.)
# PARAMETER column: "Annotation column" TYPE METACOLUMN_SEL DEFAULT description (Phenodata column to be used for the annotation)
# PARAMETER OPTIONAL use.genenames: "Represent genes with" TYPE [ids: "gene IDs", symbols: "gene names"] DEFAULT ids (Choose whether you want to print gene IDs or gene symbols in the row names of the heatmap.)
# PARAMETER OPTIONAL image.width: "Image width" TYPE INTEGER FROM 200 TO 3200 DEFAULT 1000 (Width of the image.)
# PARAMETER OPTIONAL image.height: "Image height" TYPE INTEGER FROM 200 TO 3200 DEFAULT 800 (Height of the image.)


# 12.04.2017 ML

# Functions used:
source(file.path(chipster.common.path, "tables-utils.R"))
source(file.path(chipster.common.path, "deseq-transform.R"))


## Normalise the raw counts 
# Load the counts data and extract expression value columns
dat <- read.table("data.tsv", header=T, sep="\t", row.names=1)
dat2 <- dat[,grep("chip", names(dat))]
# Call deseq_transform function:
transf.counts <- deseq_transform(dat2, method)
write.table(transf.counts, file="transformed.tsv", sep="\t", row.names=T, quote=F)

## Merge the tables 
# Read the input names to check the DE list file name, then create a new "input file table" for merging function:
input.names <- read.table("chipster-inputs.tsv", header=F, sep="\t")
listname <- as.character(input.names$V2[input.names$V1=="de-list.tsv"])
input.names2 <- data.frame(c("transformed.tsv", "de-list.tsv"), c( "transformed.tsv", listname ))
#input.names2 <- data.frame(c("data.tsv", "de-list.tsv"), c( "transformed.tsv", listname ))
colnames(input.names2) <-  c("V1", "V2")
# Call merge_tables function in tables-utils.R :
merged <- merge_tables(input.names2, "no")

# For debugging:
# write.table(list.files(), "log.txt", sep="\t", row.names=T, col.names=T, quote=F)

# filter (=select columns from) the full table for printing the table & plotting:
phenodata<-read.table("phenodata.tsv", header=T, sep="\t")
nro.samples <- nrow(phenodata)
merged_part <- merged[,c(1:nro.samples,15:20)]
dat2 <- merged[,1:nro.samples]
write.table(merged_part, "vst-transformed-counts.tsv", sep="\t", row.names=T, col.names=T, quote=F)


# Annotated heatmap

col_from<-phenodata[,pmatch(column,colnames(phenodata))] # select the column from phenodata for plotting sample names
# Use gene symbols as rownames instead of gene IDs
if (use.genenames == "symbols") {
    if (is.null(merged$symbol)){
	stop("CHIPSTER-NOTE: Sorry, there is no symbol column available in the list of differentially expressed genes. Try using the Annotation Ensembl identifiers tool in the Utilities category and use the result file as input!")
    }else{
        row_from <- merged$symbol
    }
}
pdf(file="heatmap.pdf", width=image.width/72, height=image.height/72)

# Change this to pheatmap when the package is production
heatmap(data.matrix(dat2), labCol = col_from, labRow = row_from, col = topo.colors(20)) 

# nicer version, but needs install.packages("pheatmap")
# library("pheatmap")
# pheatmap(data.matrix(dat2), labels_col = col_from, cellwidth=image.width/30 , cellheight=image.height/30)
 dev.off()

# EOF

