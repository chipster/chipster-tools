# TOOL spatial-transcriptomics-seurat-diffexp-chosen-clusters.R: "Seurat v4 -Identify spatially variable genes based on clusters" (This tool lists the differentially expressed genes between two user defined clusters.) 
# INPUT seurat_spatial_obj_pca.Robj: "Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL Markerplot.pdf
# OUTPUT OPTIONAL spatially_variable_genes.tsv
# PARAMETER OPTIONAL cluster1: "First cluster" TYPE INTEGER DEFAULT 1 (Cluster you want to identify the differentially expressed for.)
# PARAMETER OPTIONAL cluster2: "Second cluster" TYPE INTEGER DEFAULT 2 (A second cluster for comparison.)
# PARAMETER OPTIONAL min_pct: "Limit testing to genes which are expressed in at least this fraction of spots" TYPE DECIMAL DEFAULT 0.1 (Test only genes which are detected in at least this fraction of spots in either cluster. Withholding infrequently expressed genes will speed up testing.)
# PARAMETER OPTIONAL logfc_threshold: "Limit testing to genes which show at least this fold difference" TYPE DECIMAL DEFAULT 0.25 (Test only genes which show on average at least this log2 fold difference between the two groups of spots. Increasing the threshold speeds up testing, but can also miss weaker signals.)
# PARAMETER OPTIONAL test: "Test for differential expression" TYPE [wilcox: wilcox, MAST: MAST] DEFAULT wilcox
# PARAMETER OPTIONAL only_pos: "Report only positive marker genes" TYPE [FALSE, TRUE] DEFAULT FALSE (By default, this tool lists all markers. When this parameter is set to TRUE, only genes with positive log2 fold change are listed in the result file.)
# RUNTIME R-4.2.3-single-cell
# TOOLS_BIN ""

# 2022-07-29 IH
# 2022-10-20 ML Add output for spatially_variable_genes.tsv
# 2023-02-23 LG Add parameters min.pct, logfc.threshold, and only.pos

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

# Load the R-Seurat-object (called seurat_obj)
load("seurat_spatial_obj_pca.Robj")

# kokeillaan:
seurat_obj = PrepSCTFindMarkers(object = seurat_obj, assay = "SCT")


# Differential expression
de_markers <- FindMarkers(seurat_obj, ident.1 = cluster1, ident.2 = cluster2, test.use = test, logfc.threshold = logfc_threshold, min.pct = min_pct, only.pos = only_pos)

# Print out markers into a table:
# name.for.file <- paste("spatially_variable_genes_cluster", cluster1, "_vs_cluster", cluster2, ".tsv", sep="")
write.table(as.matrix(de_markers), file = "spatially_variable_genes.tsv", sep = "\t", row.names = T, col.names = T, quote = F)


# Open the pdf file for plotting
pdf(file="Markerplot.pdf", , width=9, height=12) 

# ncol = 3?
SpatialFeaturePlot(object = seurat_obj, features = rownames(de_markers)[1:3], alpha = c(0.1, 1))

dev.off() # close the pdf

#EOF


