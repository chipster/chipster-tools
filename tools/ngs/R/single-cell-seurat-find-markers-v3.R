# TOOL single-cell-seurat-find-markers-v3.R: "Seurat v3 BETA -Compare clusters" (Find the markers for a specific cluster compared to another cluster.)
# INPUT seurat_obj.Robj: "Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL log.txt
# OUTPUT OPTIONAL seurat_obj_2.Robj
# OUTPUT OPTIONAL markers.tsv
# PARAMETER OPTIONAL cluster: "Number of the cluster of interest" TYPE INTEGER DEFAULT 1 (Number of the cluster of interest.)
# PARAMETER OPTIONAL cluster2: "Cluster to compare to" TYPE STRING DEFAULT "all others" (Number\(s\) of the cluster\(s\) to compare to. By default, the cells in the cluster of interest are compared to all other cells in other clusters. You can also compare to a group of clusters, just separate the clusters with a comma \(,\) . )
# PARAMETER OPTIONAL minpct: "Min.pct" TYPE DECIMAL DEFAULT 0.25 (Only test genes that are detected in a minimum fraction of min.pct cells in either of the two populations. Meant to speed up the function by not testing genes that are very infrequently expression)
# PARAMETER OPTIONAL threshuse: "Differential expression threshold for a cluster marker gene" TYPE DECIMAL DEFAULT 0.25 (Limit testing to genes which show, on average, at least X-fold difference, in natural logarithm ln scale, between the two groups of cells. Increasing thresh.use speeds up the function, but can miss weaker signals.)
# PARAMETER OPTIONAL test.type: "Which test to use for finding marker genes" TYPE [wilcox, bimod, roc, t, tobit, poisson, negbinom, MAST, DESeq2] DEFAULT wilcox (Denotes which test to use. Seurat currently implements \"wilcox\" \(Wilcoxon rank sum test, default\), \"bimod\" \(likelihood-ratio test for single cell gene expression\), \"roc\" \(standard AUC classifier\), \"t\" \(Students t-test\), \"tobit\" \(Tobit-test for differential gene expression\), \"MAST\" \(GLM-framework that treates cellular detection rate as a covariate\), \"poisson\", and \"negbinom\". The latter two options should be used on UMI datasets only, and assume an underlying poisson or negative-binomial distribution.)
# RUNTIME R-3.6.1


# 13.06.2017 ML
# 2018-01-11 ML update Seurat version to 2.2.0
# 2020-06-17 ML update for comparing one cluster to another or more clusters


library(Seurat)
library(dplyr)
library(Matrix)
library(gplots)

# Load the R-Seurat-object (called seurat_obj)
load("seurat_obj.Robj")


# Comparing to all other cells (default): 
if (cluster2 == "all others") { 
  cluster_markers <- FindMarkers(seurat_obj, ident.1 = cluster, min.pct = minpct)
}else { 
# comparing to another user determined cluster(s):
  cluster2_fixed <- as.numeric(unlist(strsplit(cluster2, ",")))
  cluster_markers <- FindMarkers(seurat_obj, ident.1 = cluster, ident.2 = cluster2_fixed, min.pct = minpct, logfc.threshold = threshuse, test.use = test.type)
}

write.table(as.matrix(cluster_markers), file = "markers.tsv", sep = "\t", row.names = T, col.names = T, quote = F)

# Save the Robj for the next tool -not necessary here
# save(seurat_obj, file = "seurat_obj_2.Robj")

# EOF
