# TOOL spatial-transcriptomics-seurat-FindSpatiallyVariableFeatures.R: "Seurat v4 -FindSpatiallyVariableFeatures: Identification of spatially variable features" (Identify molecular features find features that have spatial patterning without taking clusters or spatial annotation into account.)
# INPUT seurat_spatial_obj_pca.Robj: "Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL Markerplot2.pdf
# PARAMETER OPTIONAL max_transparency: "Maximum transparency" TYPE DECIMAL DEFAULT 1 (Transparency of the points. Default is 1.)
# RUNTIME R-4.1.0-single-cell


library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)


# Load the R-Seurat-object (called seurat_obj)
load("seurat_spatial_obj_pca.Robj")

#Find spatially variable features using markvariogram
seurat_obj <- FindSpatiallyVariableFeatures(seurat_obj, assay = "SCT", features = VariableFeatures(seurat_obj)[1:1000], selection.method = "markvariogram")

# Open the pdf file for plotting
pdf(file="Markerplot2.pdf", , width=9, height=12) 

#Visualise the identified top features
top.features <- head(SpatiallyVariableFeatures(seurat_obj, selection.method = "markvariogram"), 6)
SpatialFeaturePlot(seurat_obj, features = top.features, ncol = 3, alpha = c(0.1, 1))

# close the pdf
dev.off()

#EOF
