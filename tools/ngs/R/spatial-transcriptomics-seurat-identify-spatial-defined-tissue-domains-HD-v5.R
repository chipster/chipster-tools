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
# RUNTIME R-4.5.1-seurat5
# SLOTS 5
# TOOLS_BIN ""

resolution = 0.5
dims.reduction = 30
label.size = 2
num_cores = as.numeric(chipster.threads.max)
parallel = TRUE

# 2026-07-15 JV


# Matrix version 1.8.0 needed for zgeMatrix. Trying to install it on the VM Runtime R-4.5.1-seurat5 through github

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(Biobase)
library("SeuratWrappers")
library("Banksy")

print("Session info:")
sessionInfo()
# BANKSY NEEDED

#remotes::install_github("prabhakarlab/Banksy@devel")
#remotes::install_github("satijalab/seurat-wrappers")


# Look into these
# k_geom param
# lambda param 

load("seurat_obj_clustering.Robj")

#Spatial.008um
DefaultAssay(seurat_obj) <- "Spatial.008um"

seurat_obj <- RunBanksy(seurat_obj, lambda = 0.8, lazy = F, verbose = T, assay = DefaultAssay(seurat_obj))

DefaultAssay(seurat_obj) <- "BANKSY"

seurat_obj <- RunPCA(seurat_obj, assay = "BANKSY", reduction.name = "pca.banksy", features = rownames(seurat_obj), npcs = 30)

seurat_obj <- FindNeighbors(seurat_obj, reduction = "pca.banksy", dims = 1:dims.reduction)

seurat_obj <- FindClusters(seurat_obj, resolution = resolution, cluster.name = "banksy_cluster")

Idents(seurat_obj) <- "banksy_cluster"

# Plotting

pdf(file = "spatiaaliplotti.pdf")

    p <- SpatialDimPlot(seurat_obj, images = "slice1.008um", group.by = "banksy_cluster", label = T, repel = T, label.size = label.size)
print(p)


    banksy_cells <- CellsByIdentities(seurat_obj)

    p1 <- SpatialDimPlot(seurat_obj, cells.highlight = banksy_cells[setdiff(names(banksy_cells), "NA")], cols.highlight = c("#FFFF00", "grey50"),
    facet.highlight = T, combine = T)

    print(p1)

dev.off()

# EOF
