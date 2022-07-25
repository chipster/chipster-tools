# TOOL spatial-transcriptomics-seurat-visualise-integration.R: "Seurat v4 -Visualise integration results" (Visualise the underlying composition of cell types in each spatial spot.)
# INPUT seurat_obj_integrated.Robj: "Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL integration_plot.pdf
# PARAMETER OPTIONAL genes: "Features to plot" TYPE STRING DEFAULT "L4" (Names of the features to plot. If you list multiple gene names, use comma as separator.)
# RUNTIME R-4.1.0-single-cell

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

# Load the R-Seurat-object (called seurat_obj)
load("seurat_obj_integrated.Robj")

genes <- trimws(unlist(strsplit(genes, ",")))

# Open the pdf file for plotting
pdf(file="integration_plot.pdf", width=13, height=7) 

# Visualise chosen features
SpatialFeaturePlot(seurat_obj, features = c(genes), pt.size.factor = 1.6, ncol = 2, crop = TRUE)

# Identify spatially variable features with the cell type prediction scores calculated in the integration
seurat_obj <- FindSpatiallyVariableFeatures(seurat_obj, assay = "predictions", selection.method = "markvariogram",
    features = rownames(seurat_obj), r.metric = 5, slot = "data")
top.clusters <- head(SpatiallyVariableFeatures(seurat_obj), 4)

# Visualise spaatially variable features
SpatialPlot(object = seurat_obj, features = top.clusters, ncol = 2)

# close the pdf
dev.off()

#EOF
