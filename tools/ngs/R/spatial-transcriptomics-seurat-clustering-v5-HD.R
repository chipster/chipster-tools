# TOOL spatial-transcriptomics-seurat-clustering-v5-HD.R: "Seurat v5 HD -Clustering" (This tool performs clustering for a single sample or multiple samples that have been combined into one Seurat object.)
# INPUT seurat_object.Robj: "Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL seurat_obj_clustering.Robj
# OUTPUT OPTIONAL clustering_plots.pdf
# PARAMETER OPTIONAL dims.reduction: "Dimensions of reduction to use" TYPE INTEGER DEFAULT 30 (Dimensions of reduction to use for clustering and UMAP.)
# PARAMETER OPTIONAL res: "Resolution for granularity for clustering" TYPE DECIMAL DEFAULT 0.8 (Resolution parameter that sets the granularity of the clustering. Higher values lead to greater number of clusters.)
# RUNTIME R-4.5.1-seurat5
# SLOTS 3
# TOOLS_BIN ""

# 2026-02 ML 

# Load seurat object (called seurat_obj)
load("seurat_object.Robj")

library(Seurat)
library(dplyr)
library(Matrix)
library(gplots)
library(ggplot2)
library(patchwork)
library(Biobase)

source(file.path(chipster.common.lib.path, "tool-utils.R"))
print(package.version("Seurat"))
documentVersion("Seurat", package.version("Seurat"))

# For developers: to be considered, when multiple sample option available:
# if ("integrated.cca" %in% Reductions(seurat_obj)) {
#     selected.reduction <- "integrated.cca"
# } else {
#     selected.reduction <- "pca"
# }

assay_names <- Assays(seurat_obj) # [1] "Spatial.008um" "Spatial.016um"

pdf(file = "clustering_plots.pdf", width = 9, height = 12)

for (i in 1:length(assay_names)) {
  
    DefaultAssay(seurat_obj) <- assay_names[i] # "Spatial.008um"
    assay_bin <- assay_names[i]
    nCount_bin <- paste("nCount_",assay_names[i], sep="")
    nFeature_bin <- paste("nFeature_",assay_names[i], sep="")
    just.bin <- sub("Spatial\\.", "", assay_names[i])
    pca_bin <- paste("pca.", just.bin, sep="")
    cluster_bin <-  paste("seurat_cluster.", just.bin, sep="") #"seurat_cluster.008um"
    umap_bin <- paste("umap.", just.bin, sep="")

    seurat_obj  <- FindNeighbors(seurat_obj, reduction = pca_bin, dims = 1:dims.reduction, verbose = FALSE)
    seurat_obj<- FindClusters(seurat_obj, resolution = res,  cluster.name = cluster_bin, verbose = FALSE,)
    seurat_obj <- RunUMAP(seurat_obj, reduction = pca_bin, reduction.name = umap_bin, dims = 1:dims.reduction, verbose = FALSE)


    # Visualization
  
    # print(DimPlot(seurat_obj, reduction = "umap", group.by = "ident"))
    # print(DimPlot(seurat_obj, reduction = "umap", group.by = "orig.ident"))
    # print(SpatialDimPlot(seurat_obj, label = TRUE, label.size = 3, ncol = 2))

    dim.plot <- DimPlot(seurat_obj, reduction = umap_bin, group.by = cluster_bin, label = TRUE, repel = T) + NoLegend()
    cluster.plot <- SpatialDimPlot(seurat_obj, group.by = cluster_bin, pt.size.factor = 1.2,  label = T, repel = T, label.size = 4) + theme(legend.position = "right")
    print(dim.plot)
    print(cluster.plot)

} # end looping for bins = assays

# Close the pdf
dev.off()

# Save the Robj for the next tool
save(seurat_obj, file = "seurat_obj_clustering.Robj")

## EOF
