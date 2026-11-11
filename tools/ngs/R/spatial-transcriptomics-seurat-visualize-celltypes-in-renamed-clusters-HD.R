# TOOL spatial-transcriptomics-seurat-visualize-celltypes-in-renamed-clusters-HD.R: "Seurat v5 HD -Visualize celltypes in renamed clusters" (Predicts cell types for Visium HD bins with RCTD using an annotated scRNA-seq reference.)
# INPUT seurat_obj_clusters_renamed.Robj: "Seurat seurat_obj" TYPE GENERIC
# OUTPUT Plots.pdf
# PARAMETER wanted_cell_types: "Cell types to plot" TYPE UNCHECKED_STRING DEFAULT "L2/3 IT" (If you list multiple cell types, please use comma\(s\) \(,\) as a separator, e.g., \"L2/3 IT\,L4\".)
# RUNTIME R-4.5.1-visium-hd
# SLOTS 4
# TOOLS_BIN ""

# 2026-07-15 JV


# Pretty much Seurat Visium HD vignette code.

# In this tool Spatial.008um assay is always used by default, for now.

library(Seurat)
library(ggplot2)
library("dplyr")
library("tidyverse")
library("Matrix")
library("SeuratWrappers")
library("spacexr")
library("ragg")
library("magick")

# Sometimes issues with this one, 3Gb should be enough.
options(future.globals.maxSize = 3000 * 1024^2) # 3000 MiB

load("seurat_obj_clusters_renamed.Robj")

# set ID to RCTD label
Idents(seurat_obj) <- 'full_first_type'

# Visualize distribution of 4 interneuron subtypes
wanted_cell_types <- wanted_cell_types
wanted_cell_types <- trimws(unlist(strsplit(wanted_cell_types, ",")))

# Keep only the cell types found in the RCTD results,
# so that CellsByIdentities() and layer_table[wanted_cell_types, ] don't fail
available_types <- unique(as.character(seurat_obj$full_first_type))
missing_types <- setdiff(wanted_cell_types, available_types)
wanted_cell_types <- intersect(wanted_cell_types, available_types)

if (length(wanted_cell_types) == 0) {
    stop(paste0(
        "None of the cell types to visualize were found in the RCTD results. ",
        "Available cell types: ", paste(sort(available_types), collapse = ", ")
    ))
}

if (length(missing_types) > 0) {
    print(paste0(
        "These cell types were not found in the RCTD results and are not plotted: ",
        paste(missing_types, collapse = ", ")
    ))
}

cell_ids <- CellsByIdentities(seurat_obj, idents = wanted_cell_types)
p <- SpatialDimPlot(seurat_obj, cells.highlight = cell_ids, cols.highlight = c("#FFFF00","grey50"), facet.highlight = T, combine=T, ncol=4)



# drop = FALSE keeps the 2D table even when only one cell type is chosen
layer_table <- table(seurat_obj$full_first_type, seurat_obj$layer_id)[wanted_cell_types, , drop = FALSE]

props <- reshape2::melt(prop.table(layer_table), margin = 1)
p2 <- ggplot(props, aes(x = Var1, y = value, fill = Var2)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(x = "Cell type", y = "Proportion", fill = "Layer") +
  theme_classic()

# Save each plot as a PNG and combine the PNGs into one PDF.
# Rasterizing keeps the PDF small and fast to open, even with hundreds of thousands of bins.
plots <- list(
    celltypes_spatial = p,
    celltype_layer_proportions = p2
)

output_dir <- file.path(getwd(), "Celltype_plots")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Number the files so that the PDF pages keep the plotting order
png_files <- file.path(output_dir, sprintf("%02d_%s.png", seq_along(plots), names(plots)))

for (i in seq_along(plots)) {
    agg_png(
        filename = png_files[i],
        width = 12,
        height = 8,
        res = 300,
        bitsize = 16,
        units = "in"
    )
    print(plots[[i]])
    dev.off()
}

images <- image_read(png_files)
image_write(images, format = "pdf", path = "Plots.pdf", density = 300)

# EOF
