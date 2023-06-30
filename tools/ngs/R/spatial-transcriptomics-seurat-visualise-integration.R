# TOOL spatial-transcriptomics-seurat-visualise-integration.R: "Seurat v4 -Visualize cell types after integration with scRNA-seq data" (Visualize the underlying composition of cell types in each spatial spot.)
# INPUT seurat_obj_integrated.Robj: "Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL integration_plot.pdf
# PARAMETER OPTIONAL genes: "Features to plot" TYPE STRING DEFAULT "L4" (Names of the features to plot. If you list multiple gene names, use comma as separator.)
# PARAMETER OPTIONAL method.to.use: "Method to use" TYPE [markvariogram: markvariogram, moransi: moransi] DEFAULT markvariogram (Method to use. Markvariogram takes longer to run, Morans I is faster.)
# PARAMETER OPTIONAL number.of.top.features: "Number of features to plot" TYPE INTEGER DEFAULT 6 (How many top features to plot.)
# RUNTIME R-4.2.3-single-cell
# TOOLS_BIN ""

# 2022-08-05 IH
# 2022-10-20 ML Add moransi option 

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

# Load the R-Seurat-object (called seurat_obj)
load("seurat_obj_integrated.Robj")

genes <- trimws(unlist(strsplit(genes, ",")))

# Identify spatially variable features with the cell type prediction scores calculated in the integration
seurat_obj <- FindSpatiallyVariableFeatures(seurat_obj, assay = "predictions", selection.method = method.to.use,
                                            features = rownames(seurat_obj), r.metric = 5, slot = "data")
top.clusters <- head(SpatiallyVariableFeatures(seurat_obj), number.of.top.features)

# Open the pdf file for plotting
pdf(file="integration_plot.pdf", width=13, height=7) 

# Visualise chosen features
SpatialFeaturePlot(seurat_obj, features = c(genes), pt.size.factor = 1.6, ncol = 2, crop = TRUE)
# Visualise spatially variable features
SpatialPlot(object = seurat_obj, features = top.clusters, ncol = 2)

# close the pdf
dev.off()

#EOF
