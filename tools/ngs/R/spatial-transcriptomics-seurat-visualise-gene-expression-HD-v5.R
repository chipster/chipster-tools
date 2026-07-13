# TOOL spatial-transcriptomics-seurat-visualise-gene-expression-HD-v5.R: "Seurat v5 -Visualize gene expression in spatial transcriptomics data" (This)
# INPUT seurat_obj_clustering.Robj: "Seurat object" TYPE GENERIC
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


## PARAMETER OPTIONAL bin: "Bin size for the spatial plot" TYPE INTEGER DEFAULT 2

# 2026-09-06 JV

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

#DefaultAssay(seurat_obj) <- sprintf("Spatial.%03dum", bin)


genes <- trimws(unlist(strsplit(genes, ",")))

rownames(seurat_obj) <- toupper(rownames(seurat_obj))

match <- intersect(genes, rownames(seurat_obj))

if (!match) {
  stop("CHIPSTER-NOTE: None of the given genes were found in the Seurat object")
}


pdf(file = "gene_expression_plot.pdf", width = width, height = height)

p1 <- SpatialFeaturePlot(seurat_obj, features = genes, keep.scale = color.scale, pt.size.factor = point.size, alpha = c(min_transparency, max_transparency))

print(p1)
print(p2)

dev.off()

# EOF


