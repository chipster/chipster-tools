# TOOL single-cell-seurat-integrate-v5.R: "Seurat v5 -Integrate multiple samples" (This tool integrates multiple samples to match shared cell types and states across dataset. The samples \/R-objects need to be named when created in the Seurat Setup tool.)
# INPUT OPTIONAL seurat_obj_combined.Robj: "Merged Seurat object to integrate" TYPE GENERIC
# OUTPUT seurat_obj_integrated.Robj
# OUTPUT OPTIONAL integrated_plot.pdf
# OUTPUT OPTIONAL aver_expr_in_clusters.tsv
# OUTPUT OPTIONAL log_normalized.tsv
# PARAMETER OPTIONAL normalisation.method: "Normalisation method used previously" TYPE [LogNormalize:"Global scaling normalization", SCT:"SCTransform"] DEFAULT LogNormalize (Which normalisation method was used in preprocessing, Global scaling normalization \(default, NormalizeData function used\) or SCTransform.)
# PARAMETER OPTIONAL anchor.identification.method: "Anchor identification method" TYPE [CCAIntegration:CCA, RPCAIntegration:RPCA] DEFAULT CCAIntegration (Which anchor identification method to use. By default, canonical correlation analysis CCA is used, but user can also decide to use the faster and more conservative reciprocal PCA approach. Check from the manual in which cases this option is recommended.)
# PARAMETER OPTIONAL CCstocompute: "Number of CCs to use in the neighbor search" TYPE INTEGER DEFAULT 30 (Which dimensions to use from the CCA to specify the neighbor search space. The neighbors are used to determine the anchors for the alignment.)
# PARAMETER OPTIONAL PCstocompute: "Number of PCs to use in the integration" TYPE INTEGER DEFAULT 30 (Number of PCs to use in the anchor weighting procedure. The anchors and their weights are used to compute the correction vectors, which allow the datasets to be integrated.)
# PARAMETER OPTIONAL num.dims: "Number of PCs to use for UMAP or TSNE" TYPE INTEGER DEFAULT 30 (Number of principal components to use.)
# PARAMETER OPTIONAL res: "Resolution for granularity" TYPE DECIMAL DEFAULT 0.8 (Resolution parameter that sets the granularity of the clustering. Increased values lead to greater number of clusters. Values between 0.6-1.2 return good results for single cell datasets of around 3K cells. For larger data sets, try higher resolution.)
# PARAMETER OPTIONAL reduction.method: "Visualisation of clusters with tSNE, UMAP or PCA" TYPE [umap:UMAP, tsne:tSNE, pca:PCA] DEFAULT umap (Which dimensionality reduction to use.)
# PARAMETER OPTIONAL point.size: "Point size in cluster plot" TYPE DECIMAL DEFAULT 0.5 (Point size for the dimensionality reduction plot.)
# PARAMETER OPTIONAL add.labels: "Add labels on top of clusters in plots" TYPE [TRUE: yes, FALSE: no] DEFAULT TRUE (Add cluster number on top of the cluster in UMAP and tSNE plots.)
# PARAMETER OPTIONAL output_aver_expr: "Give a list of average expression in each cluster" TYPE [T: yes, F: no] DEFAULT F (Returns an expression table for an 'average' single cell in each cluster.)
# RUNTIME R-4.3.2-single-cell
# SLOTS 2
# TOOLS_BIN ""

# add parameter if needed
# PARAMETER OPTIONAL ref.sample.names: "Samples to use as references" TYPE STRING DEFAULT "No references selected" (Names of the sample or samples you wish to use as references in integration, separated by comma. If you are integrating several large datasets, the tool might run out of memory. Choosing to use only some of them as references makes the integration more memory efficient and faster. Please note that the sample names here are case sensitive, so check how you typed the names of the samples when running the setup tool.)


# 2021-12-30 ML
# 2022-02-17 EK increased slots to 4
# 2022-04-19 ML increased slots to 5
# 2022-05-04 ML add RPCA option for anchor identification
# 2022-05-05 ML Rewrite the code, add option to use only part of samples as references
# 2023-02-03 ML Add 5 slots
# 2023-04-06 LG Remove 5 slots
# 2023-02-01 ML Return to the original 3 slots
# 2023-12-15 IH Update to Seurat v5

library(Seurat)
library(gplots)
library(ggplot2)
require(cowplot)
options(Seurat.object.assay.version = "v5")
options(future.globals.maxSize = 8000 * 1024^2)

# Load the R-Seurat-object
load("seurat_obj_combined.Robj")


# Compute UMAP without integration
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:CCstocompute, reduction = "pca")
seurat_obj <- FindClusters(seurat_obj, resolution = res, cluster.name = "unintegrated_clusters")
seurat_obj <- RunUMAP(seurat_obj, dims = 1:num.dims, reduction = "pca", reduction.name = "umap.unintegrated") # dims = Which dimensions to use as input features

# Perform integration:
# Not needed anymore because of merge
# When using RPCA, need to run PCA on each dataset using these features:
#if (anchor.identification.method == "rpca") {
    #seurat.objects.list <- lapply(X = seurat.objects.list, FUN = function(x) {
        #x <- ScaleData(x, features = features, verbose = FALSE)
        #x <- RunPCA(x, features = features, verbose = FALSE)
    #})
#}

if (anchor.identification.method == "CCAIntegration") {
  new.reduction = "integrated.cca"
} else if (anchor.identification.method == "RPCAIntegration") {
  new.reduction = "integrated.rpca"
}

# When using only the user listed samples as references:
# add later if needed
#if (ref.sample.names != "No references selected") {
    #ref.samples.name.list <- unlist(strsplit(ref.sample.names, ", "))
    # Go through the samples = R-objects in the list to see which ones are the reference samples.
    #ref.sample.numbers <- vector()
    #for (i in 1:length(seurat.objects.list)) {
        # Check, if the (first) sample name i:th sample is one of the names listed by user (=if there are any TRUEs)
        #if (any(seurat.objects.list[[i]]@meta.data$stim[1] == ref.samples.name.list)) {
            # if TRUE, save the i
            #ref.sample.numbers <- append(ref.sample.numbers, i)
        #}
    #}
#} else {
    #ref.sample.numbers <- NULL # if no samples are listed, NULL = all pairwise anchors are found (no reference/s)
#}

if (normalisation.method == "SCT") {
  if (length(seurat_obj@assays$SCT) > 0) {
    #seurat.objects.list <- PrepSCTIntegration(object.list = seurat_obj, anchor.features = features)
    print(anchor.identification.method)
    print(new.reduction)
    data.combined <- IntegrateLayers(object = seurat_obj, method = anchor.identification.method, normalization.method = "SCT", orig.reduction = "pca", new.reduction = new.reduction, dims = 1:PCstocompute,
      verbose = FALSE)
  } else {
    stop(paste("CHIPSTER-NOTE: ", "The data you provided hasn't been SCTransformed, please run SCTransform first or choose other parameter value."))
  }
} else {
    print(anchor.identification.method)
    print(new.reduction)
    data.combined <- IntegrateLayers(object = seurat_obj, method = anchor.identification.method, orig.reduction = "pca", new.reduction = new.reduction, dims = 1:PCstocompute,
      verbose = FALSE)
}
data.combined[["RNA"]] <- JoinLayers(data.combined[["RNA"]])
print(data.combined)

data.combined <- FindNeighbors(data.combined, reduction = new.reduction, dims = 1:CCstocompute)
data.combined <- FindClusters(data.combined, resolution = res)

# t-SNE and UMAP
# NOTE: let's do both tSNE AND UMAP so that both can be later visualized.
data.combined <- RunUMAP(data.combined, dims = 1:num.dims, reduction = new.reduction) # dims = Which dimensions to use as input features
data.combined <- RunTSNE(data.combined, dims = 1:num.dims, reduction = new.reduction) # dims = Which dimensions to use as input features

# Visualization
pdf(file = "integrated_plot.pdf", width = 13, height = 7) # open pdf
DimPlot(data.combined, reduction = "umap.unintegrated", group.by = c("stim", "seurat_clusters"), pt.size = point.size) + labs(title = "UMAP unintegrated")
DimPlot(data.combined, reduction = reduction.method, group.by = c("stim", "seurat_clusters"), pt.size = point.size, label = add.labels) + labs(title = "Chosen reduction method, integrated")
#plot_grid(p1, p2)
# Show both conditions in separate plots:
DimPlot(data.combined, reduction = reduction.method, split.by = "stim", pt.size = point.size, label = add.labels) + labs(title = "Chosen reduction method, integrated, samples")


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
save(data.combined, file = "seurat_obj_integrated.Robj")

## EOF
