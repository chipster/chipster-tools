# TOOL spatial-transcriptomics-seurat-clustering-plots-v5-HD.R: "Seurat v5 -Plot clusters" (When provided with a Visium HD Seurat object with cluster information, this tool visualizes clusters on a UMAP plot and on their spatial location. This allows you to plot clusters with new names, if you have annotated the clusters.)
# INPUT seurat_object.Robj: "Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL seurat_obj_clustering_plots.Robj
# OUTPUT OPTIONAL clustering_plots.pdf
# PARAMETER banksy: "Plot for BANKSY clusters" TYPE [TRUE: yes, FALSE: no] DEFAULT FALSE (Select yes if you identified spatially defined tissue domains with BANKSY.)
# RUNTIME R-4.5.1-visium-hd
# SLOTS 4
# TOOLS_BIN ""

# 2026-02 ML
# 2026-10 JV

# On this script in the future, one could add params for widht, height, alpha, color for highlight etc.

# Load seurat object (called seurat_obj)
load("seurat_object.Robj")

library(Seurat)
library(dplyr)
library(Matrix)
library(gplots)
library(ggplot2)
library(patchwork)
library(Biobase)
library("ragg")
library("magick")

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

output_dir <- paste0(getwd(), "/Plot_clusters")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

plots <- list()


if (banksy) {
    # Check if banksy truly exists
    if (!("BANKSY" %in% Reductions(seurat_obj) && "banksy_cluster" %in% colnames(seurat_obj@meta.data))) {
        stop("CHIPSTER-NOTE: BANKSY not found in your input data!")
    } else {
        bk_cells <- rownames(seurat_obj@meta.data)[!is.na(seurat_obj$banksy_cluster)]
        banksy_assay <- if (any(grepl("s_008um", bk_cells))) "Spatial.008um" else "Spatial.016um"

        DefaultAssay(seurat_obj) <- banksy_assay
        img_name <- Images(seurat_obj, assay = banksy_assay)
        Idents(seurat_obj) <- "banksy_cluster"

        # Not needed for not, this plot type
        #     dim.plot_banksy <- DimPlot(
        #     seurat_obj,
        #     reduction = "BANKSY",
        #     group.by = "banksy_cluster",
        #     label = TRUE,
        #     repel = T,
        #     alpha = 0.3
        # ) + ggtitle(paste0(banksy_assay," ", "BANKSY clusters"))

        p0 <- SpatialDimPlot(
            seurat_obj,
            images = img_name,
            group.by = "banksy_cluster",
            label = TRUE,
            repel = TRUE,
            label.size = 4
        ) +
            ggtitle(paste0(banksy_assay, " ", "Spatial plot of BANKSY clusters"))

        banksy_cells <- CellsByIdentities(seurat_obj, cells = bk_cells)
        banksy_cells <- banksy_cells[setdiff(names(banksy_cells), "NA")]

        p1 <- SpatialDimPlot(
            seurat_obj,
            images = img_name,
            cells.highlight = banksy_cells,
            cols.highlight = c("#FFFF00", "grey50"),
            facet.highlight = TRUE,
            combine = TRUE
        ) +
            plot_annotation(
                title = paste0(banksy_assay, " ", "BANKSY clusters highlighted on top of tissue image"),
                theme = theme(plot.title = element_text(hjust = 0.5, size = 16))
            ) &
            NoLegend()

        # If you wonder why 0_ and 1_, its because we wanna show these first
        # The plot list is in alphabetical order and numbers first

        #plots[[paste0(banksy_assay, "00_dim.plot_banksy")]] <- dim.plot_banksy
        plots[[paste0(banksy_assay, "00_p0_banksy")]] <- p0
        plots[[paste0(banksy_assay, "10_p1_banksy")]] <- p1
    }
} else {
    for (i in 1:length(assay_names)) {
        DefaultAssay(seurat_obj) <- assay_names[i] # "Spatial.008um"
        assay_bin <- assay_names[i]
        nCount_bin <- paste("nCount_", assay_names[i], sep = "")
        nFeature_bin <- paste("nFeature_", assay_names[i], sep = "")
        just.bin <- sub("Spatial\\.", "", assay_names[i])
        umap_bin <- paste0("umap.", just.bin)
        cluster_bin <- paste0("seurat_cluster.", just.bin)

        # Skip sketch assay since it is not needed after projection (Which is done in clustering tool already)
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
            alpha = 0.3
        ) +
            ggtitle(paste0(assay_bin, " ", "Projected clustering (full dataset)"))

        cluster.plot <- SpatialDimPlot(
            seurat_obj,
            group.by = cluster_bin,
            pt.size.factor = 1.2,
            label = T,
            repel = T,
            label.size = 4
        ) +
            theme(legend.position = "right") +
            ggtitle("Spatial plot of seurat clusters on top of tissue image")
        #print(dim.plot)
        #print(cluster.plot)

        # Apparently Seurat knows if an Ident string character is found in the columns so no need to use sobj@meta.data[[cluster_bin]]
        Idents(seurat_obj) <- cluster_bin
        cells <- CellsByIdentities(seurat_obj)
        cells <- cells[setdiff(names(cells), "NA")]

        highlight.plot <- SpatialDimPlot(
            seurat_obj,
            cells.highlight = cells,
            cols.highlight = c("#FFFF00", "grey50"),
            facet.highlight = TRUE,
            combine = TRUE,
            ncol = 6
        ) +
            plot_annotation(
                title = "Seurat clusters highlighted on top of tissue image",
                theme = theme(plot.title = element_text(hjust = 0.5, size = 16))
            ) &
            NoLegend()

        plots[[paste0(just.bin, "_p0")]] <- dim.plot
        plots[[paste0(just.bin, "_p1")]] <- cluster.plot
        plots[[paste0(just.bin, "_p2")]] <- highlight.plot
    }
}
# Actually plot the pngs and knit into pdf

for (i in names(plots)) {
    agg_png(
        filename = file.path(output_dir, paste0(i, ".png")),
        width = 10,
        height = 10,
        res = 300,
        bitsize = 16,
        units = "in"
    )
    print(plots[[i]])
    dev.off()
}

png_files <- list.files(
    output_dir,
    pattern = "\\.png$",
    full.names = TRUE
)

# Sort them so the PDF follows assay -> plot order
png_files <- sort(png_files)

# Read all PNGs
images <- image_read(png_files)

image_write(images, format = "pdf", path = "clustering_plots.pdf")

# EOF
