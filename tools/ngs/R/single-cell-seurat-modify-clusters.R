# TOOL single-cell-seurat-modify-clusters.R: "Seurat v4 -Remove clusters" (Remove particular clusters from data.)
# INPUT seurat_obj.Robj: "Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL log.txt
# OUTPUT seurat_obj_modified_clusters.Robj
# OUTPUT OPTIONAL ClusterPlots.pdf 
# PARAMETER cluster.identifier: "Names of the clusters to remove or include" TYPE STRING DEFAULT "3" (Name or names of the clusters you wish to keep or remove. If you list several, separate them with comma.)
# PARAMETER remove.include: "Do you wish to remove or include these clusters" TYPE [remove, include] DEFAULT remove (Select whether you want to remove or include the clusters listed above.)
# PARAMETER OPTIONAL reduction.method: "Visualisation of clusters with tSNE, UMAP or PCA" TYPE [umap:UMAP, tsne:tSNE, pca:PCA] DEFAULT umap (Which dimensionality reduction to use in plots.)
# PARAMETER OPTIONAL point.size: "Point size in cluster plot" TYPE DECIMAL DEFAULT 0.5 (Point size for the dimensionality reduction plot.)
# RUNTIME R-4.1.0-single-cell


# 2021-12-28 ML 

library(Seurat)
library(dplyr)
library(Matrix)
library(gplots)
library(ggplot2)
require(cowplot)

# Load the R-Seurat-object (called seurat_obj)
load("seurat_obj.Robj")

# QC plots, open pdf:
pdf(file="ClusterPlots.pdf",  width=13, height=7) 

if (exists("data.combined") ){
	seurat_obj <- data.combined
}

# If multiple genes are listed: (separate words from "," and remove whitespace)
if(length(grep(",", cluster.identifier)) != 0) {
   cluster.identifier <- trimws(unlist(strsplit(cluster.identifier, ",")))
}

# Sanity check: are the requested clusters available in the data:
all.clusters <- levels(seurat_obj)
match(cluster.identifier, all.clusters)
# if one of the genes is not in the list, print error message:
if (!all(!is.na(match(cluster.identifier, all.clusters)))) { 
  not.found <- cluster.identifier[is.na(match(cluster.identifier, all.clusters))==TRUE]
 #  print(paste("The gene you requested was not found in this dataset:", not.found))
  stop(paste('CHIPSTER-NOTE: ', "The cluster you listed was not found in this dataset:", not.found))
  }


if (remove.include == "include"){
    seurat_obj <- subset(x = seurat_obj, idents = cluster.identifier)
}else{ # exclude
    seurat_obj <- subset(x = seurat_obj, idents = cluster.identifier, invert = TRUE)
}

# Save the Robj for the next tool
save(seurat_obj, file="seurat_obj_modified_clusters.Robj")

# Plots to visually confirm what happened
DimPlot(seurat_obj, reduction = "tsne", label = TRUE, pt.size = point.size)
DimPlot(seurat_obj, reduction = "umap", label = TRUE, pt.size = point.size)
# Number of cells in each cluster:
cell_counts <- table(Idents(seurat_obj), seurat_obj$orig.ident)
textplot(cell_counts, halign = "center", valign = "center", cex = 1)
title(paste("Total number of cells: ", length(colnames(x = seurat_obj)), "\n Number of cells in each cluster:"))

# If there were several samples (=seurat.obj@stim is not empty) in the data, plot also sample plots:
res <- try(seurat_obj$stim, silent = TRUE)
if (class(res) != "try-error"){
  p1 <- DimPlot(seurat_obj, reduction = reduction.method, group.by = "stim", pt.size = point.size)
  p2 <- DimPlot(seurat_obj, reduction = reduction.method, pt.size = point.size, label=TRUE)
  plot_grid(p1, p2)
}
dev.off() # close the pdf

## EOF