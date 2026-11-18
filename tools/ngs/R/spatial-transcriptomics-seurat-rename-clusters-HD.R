# TOOL spatial-transcriptomics-seurat-rename-clusters-HD.R: "Seurat v5 -Rename clusters" (Renames BANKSY or Seurat clusters, for example to cell type names, and plots the renamed clusters on the tissue image.)
# INPUT seurat_obj_integrated.Robj: "Seurat seurat_obj" TYPE GENERIC
# INPUT OPTIONAL renamed_clusters.tsv: "Renamed clusters TSV" TYPE GENERIC (Optional TSV file with renamed clusters. If provided, the tool will use this file to rename the clusters instead of the default names. The tsv should contain 2 columns: cluster and new name)
# OUTPUT  seurat_obj_clusters_renamed.Robj
# OUTPUT spatial_plot_renamed_clusters.pdf
# PARAMETER clusters: "Clusters to rename" TYPE UNCHECKED_STRING DEFAULT "1,2,3,4,5" (Comma-separated list of the clusters you want to rename)
# PARAMETER names: "New names for clusters" TYPE UNCHECKED_STRING DEFAULT "Layer 2/3,Layer 4,Layer 5,Layer 6,Unknown" (Comma-separated list of new names for the selected clusters. The order of the names should match the order of the selected clusters.)
# PARAMETER OPTIONAL banksy: "Type of clusters to subset" TYPE [TRUE: BANKSY, FALSE: Seurat] DEFAULT TRUE (Seurat or BANKSY clusters. Note that seurat clusters numbers start from 0 and banksy from 1)
# RUNTIME R-4.5.1-visium-hd
# SLOTS 4
# TOOLS_BIN ""

# 2026-07-15 JV

# Pretty much Seurat Visium HD vignette code.

# In this tool Spatial.008um assay is always used by default, for now.

library(Seurat)
library(ggplot2)
library(patchwork)
library("dplyr")
library("tidyverse")
library("Matrix")
library("SeuratWrappers")
library("spacexr")
library("ragg")
library("magick")

# Sometimes issues with this one, 3Gb should be enough.
options(future.globals.maxSize = 3000 * 1024^2) # 3000 MiB

load("seurat_obj_integrated.Robj")

banksy <- as.logical(banksy)
cluster_col <- ifelse(banksy, "banksy_cluster", "seurat_cluster.projected")

if (!cluster_col %in% colnames(seurat_obj@meta.data)) {
    stop(paste0(
        "The column '",
        cluster_col,
        "' is not found in the Seurat object. ",
        ifelse(
            banksy,
            "Run BANKSY clustering first, or select no for BANKSY clusters.",
            "Run the clustering tool first."
        )
    ))
}

Idents(seurat_obj) <- cluster_col


# Make a layer_id column in the seurat_obj metadata
seurat_obj$layer_id <- 'Unknown'

# Load the tsv with new cluster names:
if (file.exists("renamed_clusters.tsv")) {
    clusternames <- read.table(
        "renamed_clusters.tsv",
        sep = "\t",
        header = TRUE,
        colClasses = "character",
        quote = "",
        comment.char = "",
        strip.white = TRUE
    )
    if (ncol(clusternames) < 2) {
        stop("renamed_clusters.tsv must have two tab-separated columns: the cluster id and the new name.")
    }
    # Use the first two columns whatever their header says
    clusternames <- clusternames[, 1:2]
    colnames(clusternames) <- c("cluster", "new_name")
} else {
    clusters <- trimws(unlist(strsplit(clusters, ",")))
    names <- trimws(unlist(strsplit(names, ",")))
    if (length(clusters) != length(names)) {
        stop(paste0(
            "The number of clusters (",
            length(clusters),
            ") and new names (",
            length(names),
            ") must be the same."
        ))
    }
    clusternames <- data.frame(cluster = clusters, new_name = names)
}

for (i in 1:nrow(clusternames)) {
    cluster <- clusternames$cluster[i]
    new_name <- clusternames$new_name[i]
    seurat_obj$layer_id[as.character(seurat_obj@meta.data[[cluster_col]]) == as.character(cluster)] <- new_name
}

# Set new Idents

Idents(seurat_obj) <- "layer_id"

p <- SpatialDimPlot(seurat_obj, group.by = "layer_id", label = TRUE)

# Highlight each renamed cluster separately
cells <- CellsByIdentities(seurat_obj)
cells <- cells[setdiff(names(cells), "NA")]

# CLAUDE: move Unknown to the end so that its panel is always the last one
cells <- cells[c(setdiff(names(cells), "Unknown"), intersect("Unknown", names(cells)))]
p1 <- SpatialDimPlot(
    seurat_obj,
    cells.highlight = cells,
    cols.highlight = c("#FFFF00", "grey50"),
    facet.highlight = TRUE,
    combine = TRUE
) +
    NoLegend() +
    plot_annotation(title = "Each renamed cluster highlighted")

# Save the plots as PNGs and combine them into a PDF, a vector PDF of HD data gets huge
output_dir <- paste0(getwd(), "/Renamed_clusters_plots")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

plots <- list(p0 = p, p1 = p1)

#plots <- rev(plots)

for (nm in names(plots)) {
    agg_png(
        filename = file.path(output_dir, paste0(nm, ".png")),
        width = 10,
        height = 10,
        res = 300,
        bitsize = 16,
        units = "in"
    )
    print(plots[[nm]])
    dev.off()
}

png_files <- sort(list.files(output_dir, pattern = "\\.png$", full.names = TRUE))

images <- image_read(png_files)
image_write(images, format = "pdf", path = "spatial_plot_renamed_clusters.pdf")

save(seurat_obj, file = "seurat_obj_clusters_renamed.Robj")

# EOF
