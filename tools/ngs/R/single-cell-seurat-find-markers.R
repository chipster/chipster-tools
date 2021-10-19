# TOOL single-cell-seurat-find-markers.R: "Seurat v4 -Find differentially expressed genes between clusters" (Find genes which are differentially expressed between given clusters or between a specified cluster and all the other cells.)
# INPUT seurat_obj.Robj: "Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL log.txt
# OUTPUT OPTIONAL seurat_obj_markers.Robj
# OUTPUT OPTIONAL markers.tsv
# PARAMETER OPTIONAL cluster: "Cluster of interest" TYPE INTEGER DEFAULT 1 (The number of the cluster of interest.)
# PARAMETER OPTIONAL cluster2: "Cluster to compare with" TYPE STRING DEFAULT "all others" (Number\(s\) of the cluster\(s\) to compare to. By default the cluster of interest is compared to cells in all other clusters. You can also compare to another cluster or a group of clusters, just separate the cluster numbers with a comma.)
# PARAMETER OPTIONAL minpct: "Limit testing to genes which are expressed in at least this fraction of cells" TYPE DECIMAL DEFAULT 0.1 (Test only genes which are detected in at least this fraction of cells in either of the two populations. Meant to speed up testing by leaving out genes that are very infrequently expressed.)
# PARAMETER OPTIONAL threshuse: "Limit testing to genes which show at least this fold difference" TYPE DECIMAL DEFAULT 0.25 (Test only genes which show on average at least this log2 fold difference, between the two groups of cells. Increasing the threshold speeds up testing, but can miss weaker signals.)
# PARAMETER OPTIONAL test.type: "Which test to use for detecting marker genes" TYPE [wilcox, DESeq2, MAST, bimod, roc, t, tobit, poisson, negbinom] DEFAULT wilcox (Seurat currently implements Wilcoxon rank sum test, bimod \(likelihood-ratio test for single cell gene expression\), roc \(standard AUC classifier\), Students t-test, Tobit-test, MAST \(GLM-framework that treates cellular detection rate as a covariate\), poisson, negbinom and DESeq2. The latter three should be used on UMI datasets only, and assume an underlying poisson or negative-binomial distribution. Note that DESeq2 is very slow and should be used only for comparisons between two clusters.)
# PARAMETER OPTIONAL only.positive: "Report only positive marker genes" TYPE [FALSE, TRUE] DEFAULT TRUE (When this parameter is set to true, only genes with positive log2 fold change are listed in the result file.)
# IMAGE comp-20.04-r-deps
# RUNTIME R-4.1.0-single-cell


# 13.06.2017 ML
# 2018-01-11 ML update Seurat version to 2.2.0
# 2020-06-17 ML update for comparing one cluster to another or more clusters
# 2020-10-12 EK add DESeq2, update min.pct default and parameter descriptions
# 2021-10-04 ML update to Seurat v4
# 2021-10-19 EK add only.pos parameter, add logfc.threshold and test.use to cluster x vs all other cells comparison

library(Seurat)
library(dplyr)
library(Matrix)
library(gplots)

# Load the R-Seurat-object (called seurat_obj)
load("seurat_obj.Robj")


# Comparing to all other cells (default): 
if (cluster2 == "all others") { 
  cluster_markers <- FindMarkers(seurat_obj, ident.1 = cluster, min.pct = minpct, logfc.threshold = threshuse, test.use = test.type, only.pos = only.positive)
}else { 
# comparing to another user determined cluster(s):
  cluster2_fixed <- as.numeric(unlist(strsplit(cluster2, ",")))
  cluster_markers <- FindMarkers(seurat_obj, ident.1 = cluster, ident.2 = cluster2_fixed, min.pct = minpct, logfc.threshold = threshuse, test.use = test.type, only.pos = only.positive)
}

write.table(as.matrix(cluster_markers), file = "markers.tsv", sep = "\t", row.names = T, col.names = T, quote = F)

# Save the Robj for the next tool -not necessary here
# save(seurat_obj, file = "seurat_obj_markers.Robj")

# EOF
