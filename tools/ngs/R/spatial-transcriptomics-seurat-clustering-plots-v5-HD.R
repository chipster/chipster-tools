# TOOL spatial-transcriptomics-seurat-clustering-plots-v5-HD.R: "Seurat v5 HD -Clustering plots" (Only plot the UMAP.)
# INPUT seurat_object.Robj: "Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL seurat_obj_clustering_plots.Robj
# OUTPUT OPTIONAL clustering_plots.pdf
# RUNTIME R-4.5.1-visium-hd
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

# mietitääns tätä sitten hieman myöhemmin:
# if ("integrated.cca" %in% Reductions(seurat_obj)) {
#     selected.reduction <- "integrated.cca"
# } else {
#     selected.reduction <- "pca"
# }
print(Assays(seurat_obj))

assay_names <- Assays(seurat_obj) # [1] "Spatial.008um" "Spatial.016um"
print(Reductions(seurat_obj))

print(colnames(seurat_obj@meta.data))

pdf(file = "clustering_plots.pdf", width = 9, height = 12)


for (i in 1:length(assay_names)) {


  
    DefaultAssay(seurat_obj) <- assay_names[i] # "Spatial.008um"
    assay_bin <- assay_names[i]
    nCount_bin <- paste("nCount_",assay_names[i], sep="")
    nFeature_bin <- paste("nFeature_",assay_names[i], sep="")
    just.bin <- sub("Spatial\\.", "", assay_names[i])
    umap_bin <- paste0("umap.", just.bin)

        if (!(umap_bin %in% Reductions(seurat_obj))) {
        next
    } 

    pca_bin <- paste("pca.", just.bin, sep="")

    if (just.bin == "sketch") {
        DefaultAssay(seurat_obj) <- "Spatial.008um"
        cluster_bin <- "seurat_cluster.projected"
        umap_bin <- "full.umap.sketch"
        pca_bin <- "full.pca.sketch"
    } else{
    cluster_bin <-  paste("seurat_cluster.", just.bin, sep="") #"seurat_cluster.008um"
    }

    # Visualization
  
    # print(DimPlot(seurat_obj, reduction = "umap", group.by = "ident"))
    # print(DimPlot(seurat_obj, reduction = "umap", group.by = "orig.ident"))
    # print(SpatialDimPlot(seurat_obj, label = TRUE, label.size = 3, ncol = 2))

print("plotting now")

    dim.plot <- DimPlot(seurat_obj, reduction = umap_bin, group.by = cluster_bin, label = TRUE, repel = T) + NoLegend()
    cluster.plot <- SpatialDimPlot(seurat_obj, group.by = cluster_bin, pt.size.factor = 1.2,  label = T, repel = T, label.size = 4) + theme(legend.position = "right")
    print(dim.plot)
    print(cluster.plot)

} 

dev.off()

#save(seurat_obj, file = "seurat_obj_clustering.Robj")

# EOF
