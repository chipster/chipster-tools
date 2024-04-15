# TOOL spatial-transcriptomics-seurat-spatially-variable-markvariogram.R: "Seurat v4 -Identify spatially variable genes without pre-annotation" (Identify genes that have spatial patterning without taking cluster information or spatial annotation into account.)
# INPUT seurat_spatial_obj_pca.Robj: "Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL Markerplot2.pdf
# OUTPUT OPTIONAL spatially_variable_genes.tsv
# PARAMETER OPTIONAL method.to.use: "Method to use" TYPE [markvariogram: markvariogram, moransi: moransi] DEFAULT markvariogram (Method to use. Mark variogram takes longer to run, Morans I is faster.)
# PARAMETER OPTIONAL number.of.top.features: "Number of features to plot" TYPE INTEGER DEFAULT 6 (How many top features to plot.)
# RUNTIME R-4.2.3-single-cell
# SLOTS 3
# TOOLS_BIN ""


# 2022-08-01 IH
# 2022-10-20 ML Add moransi option & output for spatially_variable_genes.tsv

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)


# Load the R-Seurat-object (called seurat_obj)
load("seurat_spatial_obj_pca.Robj")

# Parallelise
library(future)
plan("multisession", workers = as.integer(chipster.threads.max))

print(seurat_obj)

# Find spatially variable features using markvariogram
seurat_obj <- FindSpatiallyVariableFeatures(seurat_obj, features = VariableFeatures(seurat_obj)[1:1000], selection.method = method.to.use)

plan("default") # return to normal

# Open the pdf file for plotting
pdf(file = "Markerplot2.pdf", width = 9, height = 12)

# Visualise the identified top features
top.features <- head(SpatiallyVariableFeatures(seurat_obj, selection.method = method.to.use), number.of.top.features)

SpatialFeaturePlot(seurat_obj, features = top.features, ncol = 3, alpha = c(0.1, 1))

# close the pdf
dev.off()

# Print out markers into a table:
# name.for.file <- paste("spatially_variable_genes_cluster", cluster1, "_vs_cluster", cluster2, ".tsv", sep="")
write.table(as.matrix(top.features), file = "spatially_variable_genes.tsv", sep = "\t", row.names = T, col.names = T, quote = F)


# EOF
