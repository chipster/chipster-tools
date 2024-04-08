# TOOL spatial-transcriptomics-seurat-single-sample-analysis-v5.R: "Seurat v5 -Normalization, PCA, clustering and visualization" (This tool first normalizes data with SCTransform and detects highly variable genes. Then, the tool performs principal component analysis on the highly variable genes detected by SCTranform. It then clusters the spots and visualizes the clusters using UMAP and SpatialDimPlot. Use this tool when you want to analyze a single sample.)
# INPUT OPTIONAL seurat_object.Robj: "Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL UMAP_plot.pdf
# OUTPUT OPTIONAL seurat_spatial_obj_pca_clust.Robj
# OUTPUT OPTIONAL PCAloadings.txt
# PARAMETER OPTIONAL num.features: "Number of variable genes to return in SCTransform" TYPE INTEGER DEFAULT 3000 (Number of features to select as top variable features, i.e. how many features returned. For SCTransform, the recommended default is 3000.)
# PARAMETER OPTIONAL PCstocompute: "Number of PCs to compute" TYPE INTEGER DEFAULT 50 (Number of PCs to compute in PCA.)
# PARAMETER OPTIONAL loadings: "Print loadings in a file" TYPE [TRUE: yes, FALSE: no] DEFAULT FALSE (Print the PC loadings to a txt file.)
# PARAMETER OPTIONAL num.of.genes.loadings: "Number of genes to list in the loadings file" TYPE INTEGER DEFAULT 5 (How many genes to list in the loadings txt file.)
# PARAMETER OPTIONAL dims.reduction: "Dimensions of reduction to use as input" TYPE INTEGER DEFAULT 30 (Number of dimensions of reduction to use for clustering and UMAP.)
# PARAMETER OPTIONAL res: "Resolution for granularity for clustering" TYPE DECIMAL DEFAULT 0.8 (Resolution parameter that sets the granularity of the clustering. Increased values lead to greater number of clusters.)
# RUNTIME R-4.2.3-seurat5
# TOOLS_BIN ""

# 2022-07-27 IH
# 2022-10-20 ML Add UMAP plot with sample names
# 2024-03-21 EP Update to Seurat v5 (and add SCTransform to this tool)

library(Seurat)
library(dplyr)
library(Matrix)
library(gplots)
library(ggplot2)
library(patchwork)
library(Biobase)

print(package.version("Seurat"))
documentVersion("Seurat", package.version("Seurat"))

# Load seurat object (called seurat_obj)
load("seurat_object.Robj")

# SCTransform normalizes the data, detects high-variance features, and stores the data in the SCT assay.
seurat_obj <- SCTransform(seurat_obj, assay = "Spatial", variable.features.n = num.features, verbose = FALSE)

# PCA
seurat_obj <- RunPCA(seurat_obj, assay = "SCT", npcs = PCstocompute, verbose = FALSE)

# PCA genes in txt file
if (loadings == TRUE) {
    sink("PCAloadings.txt")
    print(seurat_obj[["pca"]], dims = 1:PCstocompute, nfeatures = num.of.genes.loadings)
    sink()
}

# Clusters
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:dims.reduction)
seurat_obj <- FindClusters(seurat_obj, resolution = res, verbose = FALSE)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:dims.reduction)

# Open the pdf file for plotting
pdf(file = "UMAP_plot.pdf", , width = 9, height = 12)

# Visualization
DimPlot(seurat_obj, reduction = "umap", group.by = "ident")
DimPlot(seurat_obj, reduction = "umap", group.by = "orig.ident")
SpatialDimPlot(seurat_obj, label = TRUE, label.size = 3)

# Close the pdf
dev.off()

# Save the Robj for the next tool
save(seurat_obj, file = "seurat_spatial_obj_pca_clust.Robj")

## EOF
