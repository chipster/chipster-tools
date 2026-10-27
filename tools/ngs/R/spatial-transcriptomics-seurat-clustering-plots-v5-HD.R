# TOOL spatial-transcriptomics-seurat-clustering-plots-v5-HD.R: "Seurat v5 -Plot clusters" (When provided with a Visium HD Seurat object with cluster information, this tool visualizes clusters on a UMAP plot and on their spatial location. This allows you to plot clusters with new names, if you have annotated the clusters.)
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
#print(Assays(seurat_obj))

#print(Reductions(seurat_obj))

#print(colnames(seurat_obj@meta.data))

# Add png to pdf method here at some point, currently pdf size can be over 0.5Gb...

assay_names <- Assays(seurat_obj) # [1] "Spatial.008um" "Spatial.016um" "sketch"


# Add also a banksy param that decides which plots are drawn -> Makes this tool more efficient

# if (banksy) {}

pdf(file = "clustering_plots.pdf", width = 9, height = 12)


for (i in 1:length(assay_names)) {
    DefaultAssay(seurat_obj) <- assay_names[i] # "Spatial.008um"
    assay_bin <- assay_names[i]
    nCount_bin <- paste("nCount_", assay_names[i], sep = "")
    nFeature_bin <- paste("nFeature_", assay_names[i], sep = "")
    just.bin <- sub("Spatial\\.", "", assay_names[i])
    umap_bin <- paste0("umap.", just.bin)
    cluster_bin <- paste0("seurat_cluster.", just.bin)
    

    if (DefaultAssay(seurat_obj) == "sketch") {
        next
    }

    if (DefaultAssay(seurat_obj) == "Spatial.008um") {
        umap_bin <- "full.umap.sketch"
        cluster_bin <- "seurat_cluster.projected"
    }

    

    dim.plot <- DimPlot(
        seurat_obj,
        reduction = umap_bin,
        group.by = cluster_bin,
        label = TRUE,
        repel = T,
        alpha = 0.1
    ) +
        NoLegend()

    cluster.plot <- SpatialDimPlot(
        seurat_obj,
        group.by = cluster_bin,
        pt.size.factor = 1.2,
        label = T,
        repel = T,
        label.size = 4
    ) +
        theme(legend.position = "right")
    print(dim.plot)
    print(cluster.plot)

    # Apparently Seurat knows if an Ident string character is found in the columns so no need to use sobj@meta.data[[cluster_bin]]
    Idents(seurat_obj) <- cluster_bin
    cells <- CellsByIdentities(seurat_obj)
    cells <- cells[setdiff(names(cells), "NA")]

    print(
        SpatialDimPlot(
            seurat_obj,
            cells.highlight = cells,
            cols.highlight = c("#FFFF00", "#aaaaaa"),
            facet.highlight = TRUE,
            combine = TRUE,
            ncol = 6
        ) +
            NoLegend()
    )
    
    if ("BANKSY" %in% Reductions(seurat_obj) &&
        "banksy_cluster" %in% colnames(seurat_obj@meta.data)) {

        bk_cells <- rownames(seurat_obj@meta.data)[!is.na(seurat_obj$banksy_cluster)]
        banksy_assay <- if (any(grepl("s_008um", bk_cells))) "Spatial.008um" else "Spatial.016um"

        DefaultAssay(seurat_obj) <- banksy_assay
        img_name <- Images(seurat_obj, assay = banksy_assay)
        Idents(seurat_obj) <- "banksy_cluster"

        p0 <- SpatialDimPlot(seurat_obj, images = img_name, group.by = "banksy_cluster",
                            label = TRUE, repel = TRUE, label.size = 4)

        banksy_cells <- CellsByIdentities(seurat_obj, cells = bk_cells)
        banksy_cells <- banksy_cells[setdiff(names(banksy_cells), "NA")]

        p1 <- SpatialDimPlot(seurat_obj, images = img_name,
                            cells.highlight = banksy_cells,
                            cols.highlight = c("#FFFF00", "grey50"),
                            facet.highlight = TRUE, combine = TRUE)
        print(p0)
        print(p1)
    }
  
}
dev.off()

# EOF
