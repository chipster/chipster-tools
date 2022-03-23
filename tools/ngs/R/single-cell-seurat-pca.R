# TOOL single-cell-seurat-pca.R: "Seurat v4 -PCA" (Principal component analysis on the highly variable genes across the single cells. The plots from this tool help you to estimate the number of principal components to be used in the clustering step.)
# INPUT OPTIONAL seurat_obj.Robj: "Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL PCAplots.pdf
# OUTPUT OPTIONAL seurat_obj_pca.Robj
# OUTPUT OPTIONAL PCAloadings.txt
# PARAMETER OPTIONAL num.of.pcas: "Number of PCs to compute" TYPE INTEGER DEFAULT 20 (How many principal components to compute and store. If you get an error message, try lowering the number. This might happen especially if you have low number of cells in your data.)
# PARAMETER OPTIONAL num.of.heatmaps: "Number of principal components to plot as heatmaps" TYPE INTEGER DEFAULT 12 (How many principal components to plot as heatmaps.)
# PARAMETER OPTIONAL loadings: "Print loadings in a file" TYPE [TRUE: yes, FALSE: no] DEFAULT FALSE (Print the PC loadings to a txt file.)
# PARAMETER OPTIONAL num.of.genes.loadings: "Number of genes to list in the loadings file" TYPE INTEGER DEFAULT 5 (How many genes to list in the loadings txt file.)
# RUNTIME R-4.1.0-single-cell


# OUTPUT OPTIONAL log.txt
# NOTE: num.of.pcas set to 20 to make runs faster, original default = 50.

# 2017-06-06 ML
# 2017-07-05 ML split into separate tool
# 2018-01-11 ML update Seurat version to 2.2.0
# 2018-09-26 ML add num.of.pcas parameter for datasets with fewer cells
# 2019-06-12 ML Seurat v3
# 2021-10-04 ML Update to Seurat v4

library(Seurat)
library(dplyr)
library(Matrix)
library(gplots)
library(ggplot2)

# Load the R-Seurat-object (called seurat_obj)
load("seurat_obj.Robj")

if (exists("data.combined") ){
	seurat_obj <- data.combined
}

# PCA
# The variable genes are used as input
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj), npcs = num.of.pcas)

# PCA genes in txt file
if (loadings == TRUE){
	sink("PCAloadings.txt")
	print(seurat_obj[["pca"]], dims = 1:num.of.pcas, nfeatures = num.of.genes.loadings)
	sink()
}

# PDF plots
pdf(file="PCAplots.pdf", , width=9, height=12) 
VizDimLoadings(seurat_obj, dims = 1:2, reduction = "pca") + ggtitle("Top 30 genes associated with PCs 1 & 2")
DimPlot(seurat_obj, reduction = "pca", group.by = "orig.ident") # orig.ident = otherwise colors based on cell cycle stages

# Need to check the number of cells at this point.
cells_left <- length(colnames(x = seurat_obj))
if (cells_left > 500) {
	DimHeatmap(seurat_obj, dims = 1, cells = 500, balanced = TRUE) #+ ggtitle("Heatmap for PC1")
	DimHeatmap(seurat_obj, dims = 1:num.of.heatmaps, cells = 500, balanced = TRUE) #+ ggtitle("Heatmaps for N first PCs")
}else{
	DimHeatmap(seurat_obj, dims = 1, cells = cells_left, balanced = TRUE) #+ ggtitle("Heatmap for PC1")
	DimHeatmap(seurat_obj, dims = 1:num.of.heatmaps, cells = cells_left, balanced = TRUE) #+ ggtitle("Heatmaps for N first PCs")
}
# fig.height=12,fig.width=9 
ElbowPlot(seurat_obj, ndims = num.of.pcas) + ggtitle("Amount of variation in the data explained by each PC")

# Number of cells:
textplot(paste("\v \v Number of \n \v \v cells: \n \v \v", length(colnames(x = seurat_obj))), halign="center", valign="center", cex=2) #, cex=0.8

dev.off() # close the pdf

# Save the Robj for the next tool
save(seurat_obj, file="seurat_obj_pca.Robj")

## EOF
