# TOOL spatial-transcriptomics-seurat-pca.R: "Seurat v4 -PCA, clustering, and visualization" (This tool performs principal component analysis on the highly variable genes detected by SCTranform. It then clusters the spots using those PCs, and visualizes the clusters using UMAP and SpatialDimPlot.)
# INPUT OPTIONAL seurat_object.Robj: "Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL UMAP_plot.pdf
# OUTPUT OPTIONAL seurat_spatial_obj_pca_clust.Robj
# OUTPUT OPTIONAL PCAloadings.txt
# PARAMETER OPTIONAL num.of.pcas: "Number of PCs to compute" TYPE INTEGER DEFAULT 30 (How many principal components to compute and use in UMAP. If you get an error message, try lowering the number. This might happen especially if you have low number of cells in your data.)
# PARAMETER OPTIONAL loadings: "Print loadings in a file" TYPE [TRUE: yes, FALSE: no] DEFAULT FALSE (Print the PC loadings to a txt file.)
# PARAMETER OPTIONAL num.of.genes.loadings: "Number of genes to list in the loadings file" TYPE INTEGER DEFAULT 5 (How many genes to list in the loadings txt file.)
# RUNTIME R-4.2.3-single-cell
# TOOLS_BIN ""

# 2022-07-27 IH
# 2022-10-20 ML Add UMAP plot with sample names

library(Seurat)
library(dplyr)
library(Matrix)
library(gplots)
library(ggplot2)
library(patchwork)

# Load seurat object (called seurat_obj)
load("seurat_object.Robj")

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
pdf(file="UMAP_plot.pdf", , width=9, height=12) 

# visualise the results of the clustering 
# Also draw UMAP plot with orig ident (handy when there are multiple samples)
# DimPlot(seurat_obj, reduction = "umap", group.by = c("ident", "orig.ident"))
# DimPlot(seurat_obj, reduction = "umap", label = TRUE)
DimPlot(seurat_obj, reduction = "umap", group.by = "ident")
DimPlot(seurat_obj, reduction = "umap", group.by = "orig.ident")

SpatialDimPlot(seurat_obj, label = TRUE, label.size = 3)

# close the pdf
dev.off()

# Save the Robj for the next tool
save(seurat_obj, file="seurat_spatial_obj_pca_clust.Robj")

## EOF
