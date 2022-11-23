# TOOL spatial-transcriptomics-seurat-diffexp-chosen-clusters.R: "Seurat v4 -Identify spatially variable genes based on clusters" (This tool lists the differentially expressed genes between two user defined clusters.) 
# INPUT seurat_spatial_obj_pca.Robj: "Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL Markerplot.pdf
# OUTPUT OPTIONAL spatially_variable_genes.tsv
# PARAMETER OPTIONAL cluster1: "First cluster" TYPE INTEGER DEFAULT 1 (Cluster you want to identify the differentially expressed for.)
# PARAMETER OPTIONAL cluster2: "Second cluster" TYPE INTEGER DEFAULT 2 (A second cluster for comparison.)
# PARAMETER OPTIONAL test: "Test for differential expression" TYPE [wilcox: wilcox, MAST: MAST] DEFAULT wilcox
# RUNTIME R-4.2.0-single-cell


# 2022-07-29 IH
# 2022-10-20 ML Add output for spatially_variable_genes.tsv

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

# Load the R-Seurat-object (called seurat_obj)
load("seurat_spatial_obj_pca.Robj")

# kokeillaan:
seurat_obj = PrepSCTFindMarkers(object = seurat_obj, assay = "SCT")


# Differential expression
de_markers <- FindMarkers(seurat_obj, ident.1 = cluster1, ident.2 = cluster2, test.use = test)

# Print out markers into a table:
# name.for.file <- paste("spatially_variable_genes_cluster", cluster1, "_vs_cluster", cluster2, ".tsv", sep="")
write.table(as.matrix(de_markers), file = "spatially_variable_genes.tsv", sep = "\t", row.names = T, col.names = T, quote = F)


# Open the pdf file for plotting
pdf(file="Markerplot.pdf", , width=9, height=12) 

# ncol = 3?
SpatialFeaturePlot(object = seurat_obj, features = rownames(de_markers)[1:3], alpha = c(0.1, 1))

dev.off() # close the pdf

#EOF


