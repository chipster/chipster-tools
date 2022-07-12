# TOOL spatial-transcriptomics-seurat-diffexp-chosen-clusters.R: "Seurat v4 -Identification of spatially variable features" (This tool lists the differentially expressed genes between two user defined clusters.) 
# INPUT seurat_spatial_obj_pca.Robj: "Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL Markerplot.pdf
# PARAMETER OPTIONAL cluster1: "First cluster" TYPE INTEGER DEFAULT 1 (Name of the cluster of which you want to identify the differentially expressed of.)
# PARAMETER OPTIONAL cluster2: "Second cluster" TYPE INTEGER DEFAULT 2 (Name of the cluster of which you want to identify the differentially expressed of.)
# PARAMETER OPTIONAL statistical.test: "Test for differential expression" TYPE [wilcox: wilcox, MAST: MAST] DEFAULT wilcox
# RUNTIME R-4.1.0-single-cell


library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)


# Load the R-Seurat-object (called seurat_obj)
load("seurat_spatial_obj_pca.Robj")

#Differential expression
de_markers <- FindMarkers(seurat_obj, ident.1 = cluster1, ident.2 = cluster2, test.use = statistical.test)

# Open the pdf file for plotting
pdf(file="Markerplot.pdf", , width=9, height=12) 

#ncol = 3?
SpatialFeaturePlot(object = seurat_obj, features = rownames(de_markers)[1:3], alpha = c(0.1, 1))

dev.off()# close the pdf

#EOF


