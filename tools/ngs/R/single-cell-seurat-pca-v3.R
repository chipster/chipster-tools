# TOOL single-cell-seurat-pca-v3.R: "Seurat v3 -PCA" (Principal component analysis on the highly variable genes across the single cells. The plots from this tool help you to estimate the number of principal components to be used in the clustering step.)
# INPUT OPTIONAL seurat_obj.Robj: "Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL PCAplots.pdf
# OUTPUT OPTIONAL seurat_obj_2.Robj
# PARAMETER OPTIONAL num.of.pcas: "Number of PCs to compute" TYPE INTEGER DEFAULT 20 (How many principal components to compute and store. If you get an error message, try lowering the number. This might happen especially if you have low number of cells in your data.)
# PARAMETER OPTIONAL num.of.heatmaps: "Number of principal components to plot as heatmaps" TYPE INTEGER DEFAULT 12 (How many principal components to plot as heatmaps.)
# RUNTIME R-3.6.1


# OUTPUT OPTIONAL log.txt
# OUTPUT OPTIONAL PCAgenes.txt
# NOTE: num.of.pcas set to 20 to make runs faster, original default = 50.

# 2017-06-06 ML
# 2017-07-05 ML split into separate tool
# 2018-01-11 ML update Seurat version to 2.2.0
# 2018-09-26 ML add num.of.pcas parameter for datasets with fewer cells
# 2019-06-12 ML Seurat v3

library(Seurat)
library(dplyr)
library(Matrix)
library(gplots)


# Load the R-Seurat-object (called seurat_obj)
load("seurat_obj.Robj")

if (exists("data.combined") ){
	seurat_obj <- data.combined
}

# PCA
# The variable genes are used as input
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj), npcs = num.of.pcas)

# PCA genes in txt file
# sink("PCAgenes.txt")
# print(seurat_obj[["pca"]], dims = 1:5, nfeatures = 5)
# sink()

# PDF plots
pdf(file="PCAplots.pdf", , width=9, height=12) 

VizDimLoadings(seurat_obj, dims = 1:2, reduction = "pca")
DimPlot(seurat_obj, reduction = "pca")

# Need to check the number of cells at this point.
cells_left <- length(colnames(x = seurat_obj))
if (cells_left > 500) {
	DimHeatmap(seurat_obj, dims = 1, cells = 500, balanced = TRUE)
	DimHeatmap(seurat_obj, dims = 1:num.of.heatmaps, cells = 500, balanced = TRUE)
}else{
	DimHeatmap(seurat_obj, dims = 1, cells = cells_left, balanced = TRUE)
	DimHeatmap(seurat_obj, dims = 1:num.of.heatmaps, cells = cells_left, balanced = TRUE)
}
# fig.height=12,fig.width=9 
ElbowPlot(seurat_obj)

# Number of cells:
textplot(paste("\v \v Number of \n \v \v cells: \n \v \v", length(colnames(x = seurat_obj))), halign="center", valign="center", cex=2) #, cex=0.8

dev.off() # close the pdf

# Save the Robj for the next tool
save(seurat_obj, file="seurat_obj_2.Robj")

## EOF
