# TOOL spatial-transcriptomics-seurat-spatially-variable-markvariogram.R: "Seurat v4 -Identify spatially variable genes without pre-annotation" (Identify genes that have spatial patterning without taking cluster information or spatial annotation into account.)
# INPUT seurat_spatial_obj_pca.Robj: "Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL Markerplot2.pdf
# OUTPUT OPTIONAL spatially_variable_genes.tsv
# PARAMETER OPTIONAL method.to.use: "Method to use" TYPE [markvariogram: markvariogram, moransi: moransi] DEFAULT markvariogram (Method to use. Mark variogram takes longer to run, Morans I is faster.)
# PARAMETER OPTIONAL number_of_var_feats: "Number of variable genes to use" TYPE INTEGER DEFAULT 1000 (Number of variable genes to use for identifying highest spatially variable genes. You can speed up the computation by choosing a smaller number of variable genes. This number should be less than or equal to the number of variable genes in the Seurat object.)
# PARAMETER OPTIONAL number.of.top.features: "Number of genes to plot" TYPE INTEGER DEFAULT 4 (How many top features to plot.)
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

# Number of variable features in Seurat object
variable_feats <- VariableFeatures(seurat_obj)

# Use number of selected variable features 
if (number_of_var_feats > length(variable_feats)) {
stop(paste("CHIPSTER-NOTE: ", "You have selected", number_of_var_feats, "to use as the number of variable features, but there are only", length(variable_feats), "in the object. Please select less than or equal to", length(variable_feats),"variable genes."))
} 
variable_feats <- variable_feats[1:number_of_var_feats]

# Find spatially variable features using markvariogram
seurat_obj <- FindSpatiallyVariableFeatures(seurat_obj, assay = "SCT", features = variable_feats, selection.method = method.to.use)

plan("default") # return to normal

# Open the pdf file for plotting
pdf(file = "Markerplot2.pdf", width = 9, height = 12)

DefaultAssay(seurat_obj) <- "SCT"

# Visualise the identified top features
top.features <- head(SpatiallyVariableFeatures(seurat_obj, selection.method = method.to.use), number.of.top.features)

SpatialFeaturePlot(seurat_obj, features = top.features, alpha = c(0.1, 1))

# close the pdf
dev.off()

# Print out markers into a table:
# name.for.file <- paste("spatially_variable_genes_cluster", cluster1, "_vs_cluster", cluster2, ".tsv", sep="")
write.table(as.matrix(top.features), file = "spatially_variable_genes.tsv", sep = "\t", row.names = T, col.names = T, quote = F)


# EOF
