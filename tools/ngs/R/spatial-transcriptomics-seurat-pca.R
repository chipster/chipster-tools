# TOOL spatial-transcriptomics-seurat-pca.R: "Seurat v4 - PCA, clustering, and visualisation" (Principal component analysis on the highly variable genes across the single cells. The plots from this tool help you to estimate the number of principal components to be used in the clustering step.)
# INPUT OPTIONAL seurat_obj_sctransform.Robj: "Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL PCAplots.pdf
# OUTPUT OPTIONAL seurat_spatial_obj_pca.Robj
# OUTPUT OPTIONAL PCAloadings.txt
# PARAMETER OPTIONAL num.of.pcas: "Number of PCs to compute" TYPE INTEGER DEFAULT 30 (How many principal components to compute and store. If you get an error message, try lowering the number. This might happen especially if you have low number of cells in your data.)
# PARAMETER OPTIONAL loadings: "Print loadings in a file" TYPE [TRUE: yes, FALSE: no] DEFAULT FALSE (Print the PC loadings to a txt file.)
# PARAMETER OPTIONAL num.of.genes.loadings: "Number of genes to list in the loadings file" TYPE INTEGER DEFAULT 5 (How many genes to list in the loadings txt file.)
# RUNTIME R-4.1.0-single-cell

library(Seurat)
library(dplyr)
library(Matrix)
library(gplots)
library(ggplot2)
library(patchwork)


# Load the R-Seurat-object (called seurat_obj)
load("seurat_obj_sctransform.Robj")

#PCA
seurat_obj <- RunPCA(seurat_obj, assay = "SCT", verbose = FALSE, npcs = num.of.pcas)

# PCA genes in txt file
if (loadings == TRUE){
	sink("PCAloadings.txt")
	print(seurat_obj[["pca"]], dims = 1:num.of.pcas, nfeatures = num.of.genes.loadings)
	sink()
}

#Clusters
seurat_obj <- FindNeighbors(seurat_obj, reduction = "pca", dims = 1:num.of.pcas)
seurat_obj <- FindClusters(seurat_obj, verbose = FALSE)
seurat_obj <- RunUMAP(seurat_obj, reduction = "pca", dims = 1:num.of.pcas)

# Open the pdf file for plotting
pdf(file="PCAplots.pdf", , width=9, height=12) 

#visualise the results of the clustering 
DimPlot(seurat_obj, reduction = "umap", label = TRUE)

SpatialDimPlot(seurat_obj, label = TRUE, label.size = 3)

# close the pdf
dev.off()

# Save the Robj for the next tool
save(seurat_obj, file="seurat_spatial_obj_pca.Robj")

## EOF