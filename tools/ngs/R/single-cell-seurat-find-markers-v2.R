# TOOL single-cell-seurat-find-markers-v2.R: "Seurat v2 -Find group specific markers" (Find the markers for a specific group.)
# INPUT seurat_obj.Robj: "Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL log.txt
# OUTPUT OPTIONAL seurat_obj_2.Robj
# PARAMETER OPTIONAL cluster: "Number of cluster of interest" TYPE INTEGER DEFAULT 1 (Number of the cluster of interest.)
# PARAMETER OPTIONAL minpct: "Min.pct" TYPE DECIMAL DEFAULT 0.25 (Only test genes that are detected in a minimum fraction of min.pct cells in either of the two populations. Meant to speed up the function by not testing genes that are very infrequently expression)
# RUNTIME R-3.4.3


# 13.06.2017 ML
# 2018-01-11 ML update Seurat version to 2.2.0


library(Seurat)
library(dplyr)
library(Matrix)
library(gplots)

# Load the R-Seurat-object (called seurat_obj)
load("seurat_obj.Robj")

#  Cluster the cells
# if biomarkers for one cluster
cluster_markers <- FindMarkers(seurat_obj, ident.1 = cluster, min.pct = minpct)


