# TOOL spatial-transcriptomics-seurat-visualise-gene-expression-HD-v5.R: "Seurat v5 -Visualize gene expression in spatial transcriptomics data" (Visualize the expression of selected genes in spatial transcriptomics data.)
# INPUT seurat_obj_clustering.Robj: "Seurat object" TYPE GENERIC (A Seurat object containing spatial transcriptomics data. This object has to be pre-processed and PCA has to be run)
# INPUT OPTIONAL genes.tsv: "Gene list in tsv format" TYPE GENERIC (A tab-separated file with a list of genes)
# OUTPUT OPTIONAL gene_expression_plot.pdf
# PARAMETER OPTIONAL genes: "Gene name\(s\)" TYPE STRING DEFAULT "Hpca, Ttr" (Name\(s\) of the gene to plot. If you list multiple gene names, use comma \(,\) as separator.)
# PARAMETER OPTIONAL point.size: "Point size in spatial feature plot" TYPE DECIMAL DEFAULT 1.6 (Point size for the plot. Default is 1.6)
# PARAMETER OPTIONAL min_transparency: "Minimum transparency" TYPE DECIMAL DEFAULT 1 (Transparency of the points. Default is 1. Transparency of points with lower expression can be downweighted with lower minimum.)
# PARAMETER OPTIONAL max_transparency: "Maximum transparency" TYPE DECIMAL DEFAULT 1 (Transparency of the points. Default is 1.)
# PARAMETER OPTIONAL color.scale: "Determine color scale based on all genes" TYPE [all:yes, feature:no] DEFAULT feature (Determine whether the color scale is based on all genes or individual genes. By default, the color scale is determined for each gene individually and may differ between genes.)
# PARAMETER OPTIONAL width: "Width of the output plot" TYPE INTEGER DEFAULT 10
# PARAMETER OPTIONAL height: "Height of the output plot" TYPE INTEGER DEFAULT 10
# RUNTIME R-4.5.1-seurat5
# SLOTS 4
# TOOLS_BIN ""

# Have to see whether we want to include the bin size here or maybe in earlier tools... Apparently some info is lacking if changed now.
## PARAMETER OPTIONAL bin: "Bin size for the spatial plot" TYPE INTEGER DEFAULT 2

# 2026-09-06 JV

#Load libraries
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(Biobase)


#packageVersion("Seurat")
#packageVersion("SeuratObject")
#print(package.version("Seurat"))
#documentVersion("Seurat", package.version("Seurat"))

# Load the R-Seurat-object (called seurat_obj)
load("seurat_obj_clustering.Robj")


# Check if the genes.tsv file exists. If so, read it and use it, else use the genes param.
if (file.exists("genes.tsv")) {
    genes <- read.table("genes.tsv", header = FALSE, sep = "\t", stringsAsFactors = FALSE)[,1]
    genes <- as.character(genes)
} else {
    genes <- unlist(strsplit(genes, ","))
}

# Make genes uppercase for easy match (User definitely knows which organism is used, so gene names can be given in any case(?))
genes <- toupper(genes)

# Remove leading and trailing whitespace from gene names
genes <- trimws(genes)

# Make rownames of seurat_obj also uppercase for easy match
rownames(seurat_obj) <- toupper(rownames(seurat_obj))

# Find common genes between seurat_obj and gene list
match <- intersect(genes, rownames(seurat_obj))

print(paste("Found", length(match), "genes in the Seurat object:"))
print(match)

# Check if any genes were found
if (length(match) == 0) {
  stop("CHIPSTER-NOTE: None of the given genes were found in the Seurat object")
}

# Plots
pdf(file = "gene_expression_plot.pdf", width = width, height = height)

p1 <- SpatialFeaturePlot(seurat_obj, features = genes, keep.scale = color.scale, pt.size.factor = point.size, alpha = c(min_transparency, max_transparency))

print(p1)

dev.off()

# EOF