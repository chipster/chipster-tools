# TOOL single-cell-seurat-integrated-analysis-v3.R: "Seurat v3 BETA -Integrated analysis of two samples" (This tool aligns the CCA subspaces and performs integrated analysis on the data. This tool can be used for two sample combined Seurat objects.) 
# INPUT OPTIONAL combined_seurat_obj.Robj: "Combined Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL integrated_plot.pdf
# OUTPUT OPTIONAL seurat_obj_combined.Robj
# PARAMETER OPTIONAL num.dims: "Number of CCs to use " TYPE INTEGER DEFAULT 20 (Number of canonical correlates to use. Use the plots from CCA tool to estimate how many you wish to use.)
# PARAMETER OPTIONAL res: "Resolution for granularity" TYPE DECIMAL DEFAULT 0.5 (Resolution parameter that sets the granularity of the clustering. Increased values lead to greater number of clusters. Values between 0.6-1.2 return good results for single cell datasets of around 3K cells. For larger data sets, try higher resolution.)
# PARAMETER OPTIONAL point.size: "Point size in tSNE plot" TYPE DECIMAL DEFAULT 0.5 (Point size for tSNE plot. )
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

# t-SNE and Clustering
 data.combined <- RunUMAP(data.combined, reduction = "pca", dims = 1:num.dims)  # dims = Which dimensions to use as input features
# data.combined <- RunTSNE(data.combined, reduction = "pca", dims = 1:num.dims) # dims = Which dimensions to use as input features
data.combined <- FindNeighbors(data.combined, reduction = "pca", dims = 1:num.dims) # dims = Dimensions of reduction to use as input
data.combined <- FindClusters(data.combined, resolution = res)

# Visualization
pdf(file="integrated_plot.pdf", , width=13, height=7)  # open pdf
 p1 <- DimPlot(data.combined, reduction = "umap", group.by = "stim", pt.size = point.size)
 p2 <- DimPlot(data.combined, reduction = "umap", label = TRUE, pt.size = point.size)
# p1 <- DimPlot(data.combined, reduction = "tsne", group.by = "stim", pt.size = point.size)
# p2 <- DimPlot(data.combined, reduction = "tsne", label = TRUE, pt.size = point.size)
plot_grid(p1, p2)
# Show both conditions in separate plots:
# DimPlot(data.combined, reduction = "tsne", split.by = "stim", pt.size = point.size)
DimPlot(data.combined, reduction = "umap", split.by = "stim", pt.size = point.size)


cell_counts <- table(Idents(data.combined), data.combined$stim)

textplot(cell_counts, halign="center", valign="center", cex=1)
title(paste("Total number of cells: ",length(colnames(x = data.combined)), "\n Number of cells in each cluster:" ) )

dev.off()


# Save the Robj for the next tool
save(data.combined, file="seurat_obj_combined.Robj")

## EOF


