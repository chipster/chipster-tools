# TOOL single-cell-seurat-clustering-v3.R: "Seurat v3 -Clustering and detection of cluster marker genes" (Clusters cells, performs tSNE and UMAP for visualization purposes, and finds positive marker genes for the clusters.)
# INPUT seurat_obj.Robj: "Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL log.txt
# OUTPUT OPTIONAL clusterPlot.pdf
# OUTPUT seurat_obj_2.Robj
# OUTPUT OPTIONAL markers.tsv
# OUTPUT OPTIONAL aver_expr_in_clusters.tsv
# PARAMETER OPTIONAL pcs_use: "Number of principal components to use" TYPE INTEGER DEFAULT 10 (How many principal components to use. User must define this based on the PCA-elbow and PCA plots from the setup tool. Seurat developers encourage to test with different parameters, and use preferably more than less PCs for downstream analysis.)
# PARAMETER OPTIONAL res: "Resolution for granularity" TYPE DECIMAL DEFAULT 0.5 (Resolution parameter that sets the granularity of the clustering. Increased values lead to greater number of clusters. Values between 0.6-1.2 return good results for single cell datasets of around 3K cells. For larger data sets, try higher resolution.)
# PARAMETER OPTIONAL perplex: "Perplexity, expected number of neighbors for tSNE plot" TYPE INTEGER DEFAULT 30 (Perplexity, expected number of neighbors. Default 30. Set to lower number if you have very few cells. Used for the tSNE visualisation of the clusters.)
# PARAMETER OPTIONAL point.size: "Point size in tSNE and UMAP plots" TYPE DECIMAL DEFAULT 1 (Point size for the cluster plots.)
# PARAMETER OPTIONAL test.type: "Which test to use for finding marker genes" TYPE [wilcox, bimod, roc, t, tobit, poisson, negbinom, MAST] DEFAULT wilcox (Denotes which test to use. Chipster currently implements wilcox \(Wilcoxon rank sum test\), bimod \(likelihood-ratio test for single cell gene expression\), roc \(standard AUC classifier\), t \(Students t-test\), tobit \(Tobit-test for differential gene expression\), MAST \(GLM-framework that treats cellular detection rate as a covariate\), poisson, and negbinom. The latter two options should be used on UMI datasets only, and assume an underlying poisson or negative-binomial distribution.)
# PARAMETER OPTIONAL minpct: "Limit testing to genes which are expressed in at least this fraction of cells" TYPE DECIMAL DEFAULT 0.1 (Test only genes which are detected in at least this fraction of cells in either of the two populations. Meant to speed up testing by leaving out genes that are very infrequently expressed.)
# PARAMETER OPTIONAL threshuse: "Limit testing to genes which show at least this fold difference" TYPE DECIMAL DEFAULT 0.25 (Test only genes which show on average at least this fold difference, in natural logarithm scale, between the two groups of cells. Increasing the threshold speeds up testing, but can miss weaker signals.)
# PARAMETER OPTIONAL output_aver_expr: "Give a list of average expression in each cluster" TYPE [T: yes, F: no] DEFAULT F (Returns an expression table for an 'average' single cell in each cluster.)
# RUNTIME R-3.6.1


# 09.06.2017 ML
# 2018-01-11 ML update Seurat version to 2.2.0
# 2018-07-25 ML removed onlypos -parameter
# 2018-09-26 ML add perplexity parameter for datasets with fewer cells
# 2019-02-18 ML add heatmap plotting & number of cells in each cluster
# 2019-04-10 ML add more test.types (wilcox, MAST, DESeq2), wilcox changed as the default
# 2019-06-15 EK change name to include marker gene detection
# 2019-06-28 EK add point size parameter for tSNE plot in the code
# 2019-06-12 ML Seurat v3
# 2019-09-09 ML UMAP
# 2019-09-23 EK add only.pos=TRUE
# 2019-11-20 EK remove DESeq2 option
# 2020-01-31 ML Add option to output average expression table
# 2020-10-11 EK update minpct to v3 default, update parameter descriptions

# for UMAP:
library(reticulate)
use_python("/opt/chipster/tools/miniconda3/envs/chipster_tools/bin/python")

library(Seurat)
library(dplyr)
library(Matrix)
library(gplots)

# Load the R-Seurat-object (called seurat_obj)
load("seurat_obj.Robj")

if (exists("data.combined")) {
  seurat_obj <- data.combined
}

# Cluster the cells

seurat_obj <- FindNeighbors(seurat_obj, dims = 1:pcs_use)
seurat_obj <- FindClusters(seurat_obj, resolution = res)


# Non-linear dimensional reduction (tSNE/UMAP) & number of cells in clusters
# NOTE: let's do both tSNE AND UMAP so that both can be later visualized.
seurat_obj <- RunTSNE(seurat_obj, dims = 1:pcs_use, do.fast = T, perplexity = perplex)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:pcs_use)

# Plot tSNE and UMAP
pdf(file = "clusterPlot.pdf")
DimPlot(seurat_obj, reduction = "umap", pt.size = point.size)
DimPlot(seurat_obj, reduction = "tsne", pt.size = point.size)
#TSNEPlot(object = seurat_obj, do.return = T, pt.size = point.size, plot.title = paste("Number of cells: ", length(colnames(x = seurat_obj))))

# Number of cells in each cluster:
cell_counts <- table(Idents(seurat_obj), seurat_obj$orig.ident)
textplot(cell_counts, halign = "center", valign = "center", cex = 1)
title(paste("Total number of cells: ", length(colnames(x = seurat_obj)), "\n Number of cells in each cluster:"))

# Find all markers
markers <- FindAllMarkers(seurat_obj, min.pct = minpct, logfc.threshold = threshuse, test.use = test.type, only.pos = TRUE)

if (length(warnings()) > 0) {
  # or !is.null(warnings())
  stop("CHIPSTER-NOTE: There was issue with FindAllMarkers functions with the selected test type, try another test!")
}

write.table(as.matrix(markers), file = "markers.tsv", sep = "\t", row.names = T, col.names = T, quote = F)

# Plot top10 genes of each cluster as a heatmap
top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

# Heatmap
DoHeatmap(object = seurat_obj, features = top10$gene) + NoLegend()
dev.off() # close the pdf

# Average expression table
# If requested, return expression for an 'average' single cell in each cluster.
if (output_aver_expr == "T") {
  aver_expr <- AverageExpression(object = seurat_obj)
  aver_expr_in_clusters <- aver_expr[[1]]
  # Write to table
  write.table(aver_expr_in_clusters, file = "aver_expr_in_clusters.tsv", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)

}

# Save the Robj for the next tool
save(seurat_obj, file = "seurat_obj_2.Robj")

# EOF
