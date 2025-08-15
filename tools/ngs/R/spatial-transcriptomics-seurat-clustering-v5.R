# TOOL spatial-transcriptomics-seurat-clustering-v5.R: "Seurat v5 -Clustering" (This tool performs clustering for a single sample or multiple samples that have been combined into one Seurat object.)
# INPUT seurat_object.Robj: "Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL seurat_obj_clustering.Robj
# OUTPUT OPTIONAL clustering_plots.pdf
# PARAMETER OPTIONAL dims.reduction: "Dimensions of reduction to use" TYPE INTEGER DEFAULT 30 (Dimensions of reduction to use for clustering and UMAP.)
# PARAMETER OPTIONAL res: "Resolution for granularity for clustering" TYPE DECIMAL DEFAULT 0.8 (Resolution parameter that sets the granularity of the clustering. Higher values lead to greater number of clusters.)
# RUNTIME R-4.2.3-seurat5
# SLOTS 3
# TOOLS_BIN ""

# 2024-05 EP Move from other tool and create own tool
# 2025-04 ML Print SpatialDimPlots in 2 columns

# Load seurat object (called seurat_obj)
load("seurat_object.Robj")

library(Seurat)
library(dplyr)
library(Matrix)
library(gplots)
library(ggplot2)
library(patchwork)
library(Biobase)

if ("integrated.cca" %in% Reductions(seurat_obj)) {
    selected.reduction <- "integrated.cca"
} else {
    selected.reduction <- "pca"
}

seurat_obj <- FindNeighbors(seurat_obj, reduction = selected.reduction, dims = 1:dims.reduction, verbose = FALSE)
seurat_obj <- FindClusters(seurat_obj, resolution = res, verbose = FALSE)
seurat_obj <- RunUMAP(seurat_obj, reduction = selected.reduction, dims = 1:dims.reduction, verbose = FALSE)

# Visualization
pdf(file = "clustering_plots.pdf", width = 9, height = 12)

print(DimPlot(seurat_obj, reduction = "umap", group.by = "ident"))
print(DimPlot(seurat_obj, reduction = "umap", group.by = "orig.ident"))
print(SpatialDimPlot(seurat_obj, label = TRUE, label.size = 3, ncol = 2))

# Close the pdf
dev.off()

# Save the Robj for the next tool
save(seurat_obj, file = "seurat_obj_clustering.Robj")

## EOF
