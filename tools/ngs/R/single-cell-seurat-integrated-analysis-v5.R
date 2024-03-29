# TOOL single-cell-seurat-integrated-analysis-v5.R: "Seurat v5 -Integrated analysis of multiple samples" (This tool performs integrated analysis on the data: clustering and visualisation of the clusters. This tool can be used for combined Seurat objects that have multiple samples in them.)
# INPUT seurat_obj_integrated.Robj: "Combined Seurat object" TYPE GENERIC
# OUTPUT seurat_obj_integrated_2.Robj
# OUTPUT OPTIONAL integrated_plot.pdf
# OUTPUT OPTIONAL aver_expr_in_clusters.tsv
# OUTPUT OPTIONAL log_normalized.tsv
# PARAMETER OPTIONAL normalisation.method: "Normalisation method used previously" TYPE [LogNormalize:"Global scaling normalization", SCT:"SCTransform"] DEFAULT LogNormalize (Which normalisation method was used in preprocessing, Global scaling normalization \(default, NormalizeData function used\) or SCTransform.)
# PARAMETER OPTIONAL num.dims: "Number of PCs to use " TYPE INTEGER DEFAULT 50 (Number of principal components to use. )
# PARAMETER OPTIONAL res: "Resolution for granularity" TYPE DECIMAL DEFAULT 0.8 (Resolution parameter that sets the granularity of the clustering. Increased values lead to greater number of clusters. Values between 0.6-1.2 return good results for single cell datasets of around 3K cells. For larger data sets, try higher resolution.)
# PARAMETER OPTIONAL reduction.method: "Visualisation of clusters with tSNE, UMAP or PCA" TYPE [umap:UMAP, tsne:tSNE, pca:PCA] DEFAULT umap (Which dimensionality reduction to use.)
# PARAMETER OPTIONAL point.size: "Point size in cluster plot" TYPE DECIMAL DEFAULT 0.5 (Point size for the dimensionality reduction plot.)
# PARAMETER OPTIONAL add.labels: "Add labels on top of clusters in plots" TYPE [TRUE: yes, FALSE: no] DEFAULT TRUE (Add cluster number on top of the cluster in UMAP and tSNE plots.)
# PARAMETER OPTIONAL output_aver_expr: "Give a list of average expression in each cluster" TYPE [T: yes, F: no] DEFAULT F (Returns an expression table for an 'average' single cell in each cluster.)
# RUNTIME R-4.3.1-single-cell
# SLOTS 2
# TOOLS_BIN ""

# To enable this option, please copy-paste this line above the #RUNTIME parameter:
# PARAMETER OPTIONAL output_norm_table: "Give a table of log-normalized values with cluster and sample information" TYPE [T: yes, F: no] DEFAULT F (Returns a table with the log-normalised UMI counts for all cells and all genes, along with the information on which sample and which cluster the cell belongs to.)


# 2018-16-05 ML
# 2019-02-18 ML add heatmap plotting & number of cells in each cluster
# 2019-04-08 ML number of cells in each cluster in each sample
# 09.07.2019 ML Seurat v3
# 2019-09-09 ML UMAP
# 2020-01-31 ML Add option to output average expression table
# 2020-02-25 ML Add option to output log-normalised values table (still commented)
# 2021-10-04 ML Update to Seurat v4
# 2922-02-21 EK Increase slots to 2 so that the average expression table is produced also with larger datasets
# 2022-07-21 ML Tune for SCTransform data
# 2023-11-29 IH remove python usage and update to Seurat v5


# UMAP uses R on default now 
#library(reticulate)
#Sys.setenv(RETICULATE_PYTHON = "/opt/chipster/tools-bin/miniconda3/envs/chipster_tools/bin/python")
# use_python("/opt/chipster/tools/miniconda3/envs/chipster_tools/bin/python")

library(Seurat)
library(gplots)
library(ggplot2)
require(cowplot)
options(Seurat.object.assay.version = "v5") 

# Load the R-Seurat-objects (called seurat_obj -that's why we need to rename them here)
load("seurat_obj_integrated.Robj")

# PCA (moved here from the combination tool)
#data.combined <- RunPCA(data.combined, npcs = num.dims, verbose = FALSE)

# t-SNE and UMAP
# NOTE: let's do both tSNE AND UMAP so that both can be later visualized.
#data.combined <- RunUMAP(data.combined, reduction = "pca", dims = 1:num.dims) # dims = Which dimensions to use as input features
#data.combined <- RunTSNE(data.combined, reduction = "pca", dims = 1:num.dims) # dims = Which dimensions to use as input features

# Clustering
# Computing nearest neighbor graph and SNN
#data.combined <- FindNeighbors(data.combined, reduction = "pca", dims = 1:num.dims) # dims = Dimensions of reduction to use as input
# Modularity Optimizer by Ludo Waltman and Nees Jan van Eck -Louvain algorithm
#data.combined <- FindClusters(data.combined, resolution = res)

# Visualization
pdf(file = "integrated_plot.pdf", width = 13, height = 7) # open pdf
#p1 <- DimPlot(data.combined, reduction = reduction.method, group.by = "stim", pt.size = point.size)
#p2 <- DimPlot(data.combined, reduction = reduction.method, pt.size = point.size, label = add.labels)
p1 <- DimPlot(data.combined  , reduction = "umap.unintegrated", group.by = c("stim", "seurat_clusters"))
p2 <- DimPlot(data.combined, reduction = "umap", group.by = c("stim", "seurat_annotations"))
p3 <- DimPlot(data.combined, reduction = "umap", split.by = "stim")
plot_grid(p1, p2, p3)

# Show both conditions in separate plots:
#DimPlot(data.combined, reduction = reduction.method, split.by = "stim", pt.size = point.size, label = add.labels)

cell_counts <- table(Idents(data.combined), data.combined$stim)
sums <- colSums(cell_counts)
cell_counts <- rbind(cell_counts, sums)

textplot(cell_counts, halign = "center", valign = "center", cex = 1)
title(paste("Total number of cells: ", length(colnames(x = data.combined)), "\n Number of cells in each cluster:"))

dev.off()

## Average expression table
## If requested, return expression for an 'average' single cell in each cluster.
# if (output_aver_expr == "T") {
#  aver_expr <- AverageExpression(object = data.combined)
#  aver_expr_in_clusters <- aver_expr[["integrated"]]
#  # Write to table
#  write.table(aver_expr_in_clusters, file = "aver_expr_in_clusters.tsv", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
# }

# Average expression table
# If requested, return expression for an 'average' single cell in each cluster.
if (output_aver_expr == "T") {
  aver_expr <- AverageExpression(object = data.combined)
  if (normalisation.method == "SCT") {
    aver_expr <- AverageExpression(object = data.combined, slot = "data", assay = "SCT")
  } else {
    aver_expr <- AverageExpression(object = data.combined)
  }

  aver_expr_in_clusters <- aver_expr[[1]]
  # aver_expr_in_clusters <- aver_expr[["integrated"]]
  # Write to table
  write.table(aver_expr_in_clusters, file = "aver_expr_in_clusters.tsv", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
}

# Normalised data + cluster + sample information table
# If requested, return a table with cells as rows and cluster + sample information and log-norm expression values for all genes in columns.
# If you want to use this option, uncomment the following section:
# if (output_norm_table == "T") {
#  norm_data <- GetAssayData(object = data.combined, slot = "data") # log-normalised "corrected" UMI counts
#  sample <- data.combined@meta.data$stim # sample information
#  cluster <- Idents(data.combined) # cluster information
#  # norm_data_table <- rbind(t(sample), t(cluster), as.matrix(norm_data)) # combine into one table
#  norm_data_table <- cbind.data.frame(sample, cluster, t(as.matrix(norm_data))) # combine into one table
#  # Write to table
#  write.table(norm_data_table, file = "log_normalized.tsv", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
# }


# Save the Robj for the next tool
save(data.combined, file = "seurat_obj_integrated_2.Robj")

## EOF