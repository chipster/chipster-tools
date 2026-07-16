# TOOL spatial-transcriptomics-seurat-identify-spatial-defined-tissue-domains-HD-v5.R: "Seurat v5 HD -Identify spatially-defined tissue domains" (This tool identifies differentially expressed genes between two user defined clusters and visualizes these genes on top of the tissue image.)
# INPUT seurat_obj_clustering.Robj: "Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL spatiaaliplotti.pdf
# OUTPUT OPTIONAL spatially_variable_genes.tsv
# PARAMETER OPTIONAL cluster1: "First cluster" TYPE INTEGER DEFAULT 1 (Cluster you want to identify the differentially expressed for.)
# PARAMETER OPTIONAL cluster2: "Second cluster" TYPE INTEGER DEFAULT 2 (A second cluster for comparison.)
# PARAMETER OPTIONAL min_pct: "Limit testing to genes which are expressed in at least this fraction of spots" TYPE DECIMAL DEFAULT 0.01 (Test only genes which are detected in at least this fraction of spots in either cluster. Withholding infrequently expressed genes will speed up testing.)
# PARAMETER OPTIONAL logfc_threshold: "Limit testing to genes which show at least this fold difference" TYPE DECIMAL DEFAULT 0.1 (Test only genes which show on average at least this log2 fold difference between the two groups of spots. Increasing the threshold speeds up testing, but can also miss weaker signals.)
# PARAMETER OPTIONAL test: "Test for differential expression" TYPE [wilcox: wilcox, MAST: MAST] DEFAULT wilcox
# PARAMETER OPTIONAL only_pos: "Report only positive marker genes" TYPE [FALSE, TRUE] DEFAULT FALSE (By default, this tool lists all markers. When this parameter is set to TRUE, only genes with positive log2 fold change are listed in the result file.)
# PARAMETER OPTIONAL no_of_feats: "Number of spatially variable genes to visualize" TYPE INTEGER DEFAULT 3 (Choose the number of highest variable genes to visualize.)
# PARAMETER OPTIONAL color.scale: "Determine color scale based on all genes" TYPE [all:yes, feature:no] DEFAULT feature (Determine whether the color scale is based on all genes or individual genes. By default, the color scale is determined for each gene individually and may differ between genes.)
# PARAMETER OPTIONAL label.size: "determine the label size of the plots" TYPE INTEGER DEFAULT 3
# PARAMETER chosen_clusters: "Choose clusters to be subsetted" TYPE INTEGER DEFAULT 0
# RUNTIME R-4.5.1-seurat5
# SLOTS 5
# TOOLS_BIN ""

library(Seurat)
library(SeuratObject)

load("seurat_obj_clustering.R")

Idents(seurat_obj) <- seurat_obj$seurat_clusters

chosen_clusters # A parameter

subset_obj <- subset(seurat_obj, idents = chosen_clusters)


pdf(file = "Spatial_plot.pdf")
p <- SpatialDimPlot(subset_obj, label = T, crop = T, label.size = 3) + NoLegend() + 
  SpatialDimPlot(subset_obj, label = T, crop = F, label.size = 3)

print(p)

dev.off()