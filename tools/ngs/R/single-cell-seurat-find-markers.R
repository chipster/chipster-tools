# TOOL single-cell-seurat-find-markers.R: "Seurat v4 -Find differentially expressed genes between clusters" (Find genes which are differentially expressed between given clusters or between a specified cluster and all the other cells.)
# INPUT seurat_obj.Robj: "Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL seurat_obj_markers.Robj
# OUTPUT OPTIONAL log.txt
# OUTPUT OPTIONAL markers.tsv
# OUTPUT OPTIONAL all_markers.tsv
# OUTPUT OPTIONAL Top10Heatmap.pdf
# PARAMETER OPTIONAL find.all.markers: "Find all markers" TYPE [FALSE, TRUE] DEFAULT FALSE (Give as an output a large table with markers for all the clusters. Each cluster is compared to all the other clusters. This parameter overwrites the two cluster number parameters below. You will want to filter this table with the tool in Utilities category.)
# PARAMETER OPTIONAL cluster: "Cluster of interest" TYPE INTEGER DEFAULT 1 (The number of the cluster of interest.)
# PARAMETER OPTIONAL cluster2: "Cluster to compare with" TYPE STRING DEFAULT "all others" (Number\(s\) of the cluster\(s\) to compare to. By default the cluster of interest is compared to cells in all other clusters. You can also compare to another cluster or a group of clusters, just separate the cluster numbers with a comma.)
# PARAMETER OPTIONAL minpct: "Limit testing to genes which are expressed in at least this fraction of cells" TYPE DECIMAL DEFAULT 0.1 (Test only genes which are detected in at least this fraction of cells in either of the two populations. Meant to speed up testing by leaving out genes that are very infrequently expressed.)
# PARAMETER OPTIONAL threshuse: "Limit testing to genes which show at least this fold difference" TYPE DECIMAL DEFAULT 0.25 (Test only genes which show on average at least this log2 fold difference, between the two groups of cells. Increasing the threshold speeds up testing, but can miss weaker signals.)
# PARAMETER OPTIONAL test.type: "Which test to use for detecting marker genes" TYPE [wilcox, DESeq2, bimod, roc, t, tobit, poisson, negbinom] DEFAULT wilcox (Seurat currently implements Wilcoxon rank sum test, bimod \(likelihood-ratio test for single cell gene expression\), roc \(standard AUC classifier\), Students t-test, Tobit-test, poisson, negbinom and DESeq2. The latter three should be used on UMI datasets only, and assume an underlying poisson or negative-binomial distribution. Note that DESeq2 is very slow and should be used only for comparisons between two clusters.)
# PARAMETER OPTIONAL only.positive: "Report only positive marker genes" TYPE [FALSE, TRUE] DEFAULT TRUE (When this parameter is set to true, only genes with positive log2 fold change are listed in the result file. NOTE, for listing all markers, this is currently set to FALSE regardless what you choose here.)
# PARAMETER OPTIONAL returnthresh: "p-value threshold" TYPE DECIMAL DEFAULT 0.01 (Only return markers that have a p-value < return.thresh, or a power > return.thresh, if the test is ROC)
# RUNTIME R-4.2.3-single-cell
# TOOLS_BIN ""

# 13.06.2017 ML
# 2018-01-11 ML update Seurat version to 2.2.0
# 2020-06-17 ML update for comparing one cluster to another or more clusters
# 2020-10-12 EK add DESeq2, update min.pct default and parameter descriptions
# 2021-10-04 ML update to Seurat v4
# 2021-10-19 EK add only.pos parameter, add logfc.threshold and test.use to cluster x vs all other cells comparison
# 2022-12-28 ML Add returnthresh as parameter (default = 0.01, so didn't work if you wanted all genes listed). onlypos won't work for FindAllMarkers, for some reason.
# 2023-02-01 LG Removed MAST option.
# PARAMETER OPTIONAL test.type: "Which test to use for detecting marker genes" TYPE [wilcox, DESeq2, MAST, bimod, roc, t, tobit, poisson, negbinom] DEFAULT wilcox (Seurat currently implements Wilcoxon rank sum test, bimod \(likelihood-ratio test for single cell gene expression\), roc \(standard AUC classifier\), Students t-test, Tobit-test, MAST \(GLM-framework that treates cellular detection rate as a covariate\), poisson, negbinom and DESeq2. The latter three should be used on UMI datasets only, and assume an underlying poisson or negative-binomial distribution. Note that DESeq2 is very slow and should be used only for comparisons between two clusters.)
# 2023-02-03 ML Add 5 slots
# 2023-04-06 LG Remove 5 slots

# PARAMETER OPTIONAL mincellsfeat: "Minimum number of cells expressing the feature" TYPE INTEGER DEFAULT 3 (Minimum number of cells expressing the feature in at least one of the two groups, currently only used for poisson and negative binomial tests.)
# PARAMETER OPTIONAL mincellsgroup: "Minimum number of cells" TYPE INTEGER DEFAULT 3 (Minimum number of cells in one of the groups.)


library(Seurat)
library(dplyr)
library(Matrix)
library(gplots)
library(ggplot2)

# Load the R-Seurat-object (called seurat_obj)
load("seurat_obj.Robj")

if (exists("data.combined")) {
  seurat_obj <- data.combined
}

# If FindAllMarkers:
if (find.all.markers == TRUE) {
  markers <- FindAllMarkers(seurat_obj, min.pct = minpct, logfc.threshold = threshuse, test.use = test.type, return.thresh = returnthresh, only.pos = as.logical(only.positive)) # min.cells.feature = mincellsfeat, min.cells.group = mincellsgroup,

  if (length(warnings()) > 0) {
    # or !is.null(warnings())
    stop("CHIPSTER-NOTE: There was issue with FindAllMarkers functions with the selected test type, try another test!")
  }

  # Filter based on adj-p-val (return.thresh parameter looks at the p_val column, not p_val_adj):
  markers_filtered <- markers[markers$p_val_adj<returnthresh, ]

  write.table(as.matrix(markers_filtered), file = "all_markers.tsv", sep = "\t", row.names = T, col.names = T, quote = F)

  # Plot top10 genes of each cluster as a heatmap
  top10 <- markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC)

  plot_heatmap <- DoHeatmap(object = seurat_obj, features = top10$gene, angle = 0, size = 2, hjust = 0.5) #+ NoLegend()
  ggplot2::ggsave(filename = "Top10Heatmap.pdf", plot = plot_heatmap)
} else { # Comparing only one cluster to all others or
  # Comparing to all other cells (default):
  if (cluster2 == "all others") {
    cluster_markers <- FindMarkers(seurat_obj, ident.1 = cluster, min.pct = minpct, logfc.threshold = threshuse, test.use = test.type, only.pos = only.positive, return.thresh = returnthresh)
  } else {
    # comparing to another user determined cluster(s):
    cluster2_fixed <- as.numeric(unlist(strsplit(cluster2, ",")))
    cluster_markers <- FindMarkers(seurat_obj, ident.1 = cluster, ident.2 = cluster2_fixed, min.pct = minpct, logfc.threshold = threshuse, test.use = test.type, only.pos = only.positive, return.thresh = returnthresh)
  }
  # Filter based on adj-p-val (no return.thresh parameter):
    cluster_markers_filtered <- cluster_markers[cluster_markers$p_val_adj<returnthresh, ]

  write.table(as.matrix(cluster_markers_filtered), file = "markers.tsv", sep = "\t", row.names = T, col.names = T, quote = F)

}
# Save the Robj for the next tool -not necessary here
# save(seurat_obj, file = "seurat_obj_markers.Robj")

# EOF