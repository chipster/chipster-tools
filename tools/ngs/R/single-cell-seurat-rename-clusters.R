# TOOL single-cell-seurat-rename-clusters.R: "Seurat v4 -Rename clusters" (You can use this tool to rename the clusters.)
# INPUT seurat_obj.Robj: "Seurat object" TYPE GENERIC
# INPUT cluster_names.tsv: "Cluster name table in tsv format" TYPE GENERIC
# OUTPUT OPTIONAL clusterPlotRenamed.pdf
# OUTPUT OPTIONAL seurat_obj_renamed.Robj
# PARAMETER OPTIONAL point.size: "Point size in tSNE and UMAP plots" TYPE DECIMAL DEFAULT 1 (Point size for the cluster plots.)
# RUNTIME R-4.1.0-single-cell


# OUTPUT OPTIONAL log.txt

# 2021-12-27 ML

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

# Load the tsv with new cluster names:
clusternames <- read.table("cluster_names.tsv", sep="\t", header=T, row.names=1)

# Get the names from the second column:
new.cluster.ids <- clusternames[,1]

# Some checks:
if(length(new.cluster.ids) != length(levels(seurat_obj)) ) {
    stop("CHIPSTER-NOTE: You need to give as input as many cluster names as there are clusters.")
}
if(row.names(clusternames) != levels(seurat_obj) ){
    stop("CHIPSTER-NOTE: The cluster numbers in the input table need to be the same as in the Seurat object.")

}

 #new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",  "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(seurat_obj)
seurat_obj <- RenameIdents(seurat_obj, new.cluster.ids)

# Plot with renamed clusters:
pdf(file = "clusterPlotRenamed.pdf")
DimPlot(seurat_obj, reduction = "umap", label = TRUE, pt.size = point.size) + NoLegend()
DimPlot(seurat_obj, reduction = "tsne", label = TRUE, pt.size = point.size)
# Number of cells in each cluster:
cell_counts <- table(Idents(seurat_obj), seurat_obj$orig.ident)
textplot(cell_counts, halign = "center", valign = "center", cex = 1)
title(paste("Total number of cells: ", length(colnames(x = seurat_obj)), "\n Number of cells in each cluster:"))

dev.off() # close the pdf

# Save the Robj for the next tool
save(seurat_obj, file = "seurat_obj_renamed.Robj")

# EOF