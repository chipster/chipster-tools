# TOOL single-cell-seurat-integrated-analysis-v3.R: "Seurat v3 -Integrated analysis of two samples" (This tool performs integrated analysis on the data: clustering and visualisation of the clusters. This tool can be used for two sample combined Seurat objects.) 
# INPUT OPTIONAL combined_seurat_obj.Robj: "Combined Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL integrated_plot.pdf
# OUTPUT OPTIONAL seurat_obj_combined.Robj
# PARAMETER OPTIONAL num.dims: "Number of PCs to use " TYPE INTEGER DEFAULT 20 (Number of principal components to use. )
# PARAMETER OPTIONAL res: "Resolution for granularity" TYPE DECIMAL DEFAULT 0.5 (Resolution parameter that sets the granularity of the clustering. Increased values lead to greater number of clusters. Values between 0.6-1.2 return good results for single cell datasets of around 3K cells. For larger data sets, try higher resolution.)
# PARAMETER OPTIONAL point.size: "Point size in tSNE plot" TYPE DECIMAL DEFAULT 0.5 (Point size for dimensionality plot. )
# PARAMETER OPTIONAL reduction.method: "Visualisation with tSNE, UMAP or PCA" TYPE [umap:UMAP, tsne:tSNE, pca:PCA] DEFAULT umap (Which dimensionality reduction to use.)
# RUNTIME R-3.6.1



# 2018-16-05 ML
# 2019-02-18 ML add heatmap plotting & number of cells in each cluster
# 2019-04-08 ML number of cells in each cluster in each sample
# 09.07.2019 ML Seurat v3
# 2019-09-09 ML UMAP


# for UMAP:
library(reticulate)
use_python("/opt/chipster/tools/miniconda3/envs/chipster_tools/bin/python")

library(Seurat)
library(gplots)
library(ggplot2)
require(cowplot)


# Load the R-Seurat-objects (called seurat_obj -that's why we need to rename them here)
load("combined_seurat_obj.Robj")
# combined_seurat_obj <- data.combined

# t-SNE and UMAP 
# NOTE: let's do both tSNE AND UMAP so that both can be later visualized.
data.combined <- RunUMAP(data.combined, reduction = "pca", dims = 1:num.dims)  # dims = Which dimensions to use as input features
data.combined <- RunTSNE(data.combined, reduction = "pca", dims = 1:num.dims) # dims = Which dimensions to use as input features

# Clustering
# Computing nearest neighbor graph and SNN
data.combined <- FindNeighbors(data.combined, reduction = "pca", dims = 1:num.dims) # dims = Dimensions of reduction to use as input
# Modularity Optimizer by Ludo Waltman and Nees Jan van Eck -Louvain algorithm
data.combined <- FindClusters(data.combined, resolution = res)

# Visualization
pdf(file="integrated_plot.pdf", , width=13, height=7)  # open pdf
 p1 <- DimPlot(data.combined, reduction = reduction.method, group.by = "stim", pt.size = point.size)
 p2 <- DimPlot(data.combined, reduction = reduction.method, label = TRUE, pt.size = point.size)
plot_grid(p1, p2)
# Show both conditions in separate plots:
DimPlot(data.combined, reduction = reduction.method, split.by = "stim", pt.size = point.size)

cell_counts <- table(Idents(data.combined), data.combined$stim)

textplot(cell_counts, halign="center", valign="center", cex=1)
title(paste("Total number of cells: ",length(colnames(x = data.combined)), "\n Number of cells in each cluster:" ) )

dev.off()


# Save the Robj for the next tool
save(data.combined, file="seurat_obj_combined.Robj")

## EOF



