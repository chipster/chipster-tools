# TOOL spatial-transcriptomics-seurat-integration-with-scRNA-HD-v5.R: "Seurat v5 HD -Integration with scRNA-seq data" (Predicts cell types for Visium HD bins with RCTD using an annotated scRNA-seq reference.)
# INPUT seurat_obj_clustering.Robj: "Seurat object" TYPE GENERIC
# INPUT scRNAseq_ref.Rds: "Seurat scRNA data" TYPE GENERIC
# OUTPUT RCTD_cell_type_plots.pdf
# OUTPUT OPTIONAL seurat_obj_integrated.Robj
# PARAMETER known_types: "Cell types to highlight in the plots" TYPE UNCHECKED_STRING DEFAULT "Astro,Oligo,Endo"
# PARAMETER OPTIONAL celltype_colnames: "Column name of the cell type information in the scRNA-seq Seurat object" TYPE STRING DEFAULT "subclass_label"
# PARAMETER OPTIONAL dims.reduction: "Number of dimensions to use for PCA and UMAP" TYPE INTEGER FROM 1 TO 100 DEFAULT 30 ()
# RUNTIME R-4.5.1-visium-hd
# SLOTS 10
# TOOLS_BIN ""

# 2026-07-15 JV


# Pretty much Seurat Visium HD vignette code.

# In this tool Spatial.008um assay is always used by default, for now.
assay <- as.character("Spatial.008um")
dims.reduction <- dims.reduction
num_cores <- as.numeric(chipster.threads.max)

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

load("seurat_obj_clustering.Robj")

print("Loading reference object")
ref <- readRDS("scRNAseq_ref.Rds")

print("Sketching the data")

# The input already has a "sketch" assay from the clustering tool, possibly made before
# the data was subset. Make a new sketch under new names, so that ProjectData() doesn't
# reuse the old full.pca.sketch, umap.sketch or the neighbors it has cached in the object.
# The clustering results (seurat_cluster.projected, full.umap.sketch) stay untouched.
sketch.ncells <- 50000
sketch_assay <- "sketch.rctd"
pca_sketch <- "pca.rctd"
umap_sketch <- "umap.rctd"
full_pca <- "full.pca.rctd" # created by ProjectData()
# ProjectData() names the projected UMAP paste0("full.", umap.model) -> "full.umap.rctd"
dims <- 1:dims.reduction

# GetTissueCoordinates() below takes the image of the default assay
DefaultAssay(seurat_obj) <- assay

# A subset may have fewer bins than the target. Then use all of them.
n_bins <- ncol(seurat_obj[[assay]])
ncells_use <- min(sketch.ncells, n_bins)
print(paste0("Sketching ", ncells_use, " of ", n_bins, " bins"))

# Same steps as the clustering tool: LeverageScore needs the scaled data,
# and it must be computed on the current variable features
seurat_obj <- FindVariableFeatures(seurat_obj, verbose = FALSE)
seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
seurat_obj <- SketchData(
    object = seurat_obj,
    assay = assay,
    ncells = ncells_use,
    method = "LeverageScore",
    sketched.assay = sketch_assay,
    features = VariableFeatures(seurat_obj)
)

# PCA and UMAP model of the sketched bins only
DefaultAssay(seurat_obj) <- sketch_assay
seurat_obj <- FindVariableFeatures(seurat_obj, verbose = FALSE)
seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
seurat_obj <- RunPCA(seurat_obj, assay = sketch_assay, reduction.name = pca_sketch, verbose = FALSE)
# return.model = TRUE so that ProjectData() can project the UMAP to all bins
seurat_obj <- RunUMAP(
    seurat_obj,
    reduction = pca_sketch,
    reduction.name = umap_sketch,
    return.model = TRUE,
    dims = dims,
    verbose = FALSE
)

gc()

# The exact column name of the cell types in the reference, e.g. "celltype"
celltype_colnames <- as.character(celltype_colnames)
print(celltype_colnames)

if (celltype_colnames %in% colnames(ref@meta.data)) {
    Idents(ref) <- celltype_colnames
} else {
    stop(paste0(
        "The celltype_colname '",
        celltype_colnames,
        "' is not found in the scRNA-seq Seurat object. Please check the column names of the object and provide a valid column name for cell types."
    ))
}

print(paste0("Idents of the ref: ", head(Idents(ref), 10)))

counts <- ref[["RNA"]]$counts
cluster <- as.factor(ref[[celltype_colnames, drop = TRUE]])
nUMI <- ref$nCount_RNA
levels(cluster) <- gsub("/", "-", levels(cluster))
cluster <- droplevels(cluster)

# create the RCTD reference object
reference <- Reference(counts, cluster, nUMI)

# Free up memory
rm(ref, counts, cluster, nUMI)
gc()

# Run RCTD on the sketched bins
counts_hd <- seurat_obj[[sketch_assay]]$counts
seurat_obj_cells_hd <- colnames(seurat_obj[[sketch_assay]])
coords_raw <- GetTissueCoordinates(seurat_obj)[seurat_obj_cells_hd, 1:2]

query <- SpatialRNA(coords_raw, counts_hd, colSums(counts_hd))

print("Running RCTD")
# run RCTD with min of 25 cells per cell type
RCTD <- create.RCTD(query, reference, max_cores = num_cores, CELL_MIN_INSTANCE = 25, test_mode = FALSE)
RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
# add results back to Seurat object
seurat_obj <- AddMetaData(seurat_obj, metadata = RCTD@results$results_df)

rm(RCTD, query, reference, counts_hd, coords_raw)
gc()

# Bins without an RCTD result are labelled "Unknown"
seurat_obj$first_type <- as.character(seurat_obj$first_type)
seurat_obj$first_type[is.na(seurat_obj$first_type)] <- 'Unknown'

# Project the RCTD labels from the sketched bins to all bins. full.pca.rctd doesn't exist yet,
# so all bins are projected to the new PCA. recompute.* makes sure that neighbors cached
# by an earlier ProjectData() call are never reused.
seurat_obj <- ProjectData(
    object = seurat_obj,
    assay = assay,
    sketched.assay = sketch_assay,
    sketched.reduction = pca_sketch,
    full.reduction = full_pca,
    umap.model = umap_sketch,
    dims = dims,
    refdata = list(full_first_type = "first_type"),
    recompute.neighbors = TRUE,
    recompute.weights = TRUE
)

DefaultAssay(seurat_obj) <- assay
Idents(seurat_obj) <- "full_first_type"

print("Plotting")

# Collect the plots in a list, save each as a PNG and combine the PNGs into one PDF.
# Rasterizing keeps the PDF small and fast to open, even with hundreds of thousands of bins.
plots <- list()

p <- SpatialDimPlot(seurat_obj)
plots[["spatial"]] <- p

cells <- CellsByIdentities(seurat_obj)

known_types <- known_types
known_types <- trimws(unlist(strsplit(known_types, ",")))

missing_types <- setdiff(known_types, names(cells))
known_types <- intersect(known_types, names(cells))

if (length(known_types) == 0) {
    (paste0(
        "None of the cell types to highlight were found in the RCTD results. ",
        "Available cell types: ", paste(sort(names(cells)), collapse = ", ")
    ))
}

if (length(missing_types) > 0) {
    print(paste0(
        "These cell types were not found in the RCTD results and are not highlighted: ",
        paste(missing_types, collapse = ", ")
    ))
}
print(known_types)

p_highlight3 <- SpatialDimPlot(
    seurat_obj,
    cells.highlight = cells[known_types],
    cols.highlight = c("#FFFF00", "black"),
    facet.highlight = TRUE,
    combine = TRUE,
    ncol = 4,
    image.alpha = 1
)
plots[["highlight_known_types_image"]] <- p_highlight3


# In this maybe "Unknown" could also be "NA"
p_highlight2 <- SpatialDimPlot(
    seurat_obj,
    cells.highlight = cells[setdiff(names(cells), "Unknown")],
    cols.highlight = c("#FFFF00", "grey50"),
    facet.highlight = T,
    combine = T,
    ncol = 8
)

plots[["highlight_all_types"]] <- p_highlight2

# Distribution of the cell types across the BANKSY clusters (Seurat vignette code),
# or across the Seurat clusters if BANKSY hasn't been run
if ("BANKSY" %in% Reductions(seurat_obj) && "banksy_cluster" %in% colnames(seurat_obj@meta.data)) {
    plot_cell_types_banksy <- function(data, label) {
        p <- ggplot(data, aes(x = get(label), y = n, fill = full_first_type)) +
            geom_bar(stat = "identity", position = "stack") +
            geom_text(
                aes(label = ifelse(n >= min_count_to_show_label, full_first_type, "")),
                position = position_stack(vjust = 0.5),
                size = 2
            ) +
            xlab(label) +
            ylab("# of Spots") +
            ggtitle(paste0("Distribution of Cell Types across ", label)) +
            theme_minimal()
    }

    cell_type_banksy_counts <- seurat_obj[[]] %>%
        dplyr::filter(full_first_type %in% known_types) %>%
        dplyr::count(full_first_type, banksy_cluster)

    min_count_to_show_label <- 20

    p_banksy <- plot_cell_types_banksy(cell_type_banksy_counts, "banksy_cluster")
    plots[["banksy_cell_types"]] <- p_banksy
} else {
    plot_cell_types_seurat <- function(data, label) {
        p <- ggplot(data, aes(x = get(label), y = n, fill = full_first_type)) +
            geom_bar(stat = "identity", position = "stack") +
            geom_text(
                aes(label = ifelse(n >= min_count_to_show_label, full_first_type, "")),
                position = position_stack(vjust = 0.5),
                size = 2
            ) +
            xlab(label) +
            ylab("# of Spots") +
            ggtitle(paste0("Distribution of Cell Types across ", label)) +
            theme_minimal()
    }

    cell_type_seurat_counts <- seurat_obj[[]] %>%
        dplyr::filter(full_first_type %in% known_types) %>%
        dplyr::count(full_first_type, seurat_cluster.projected)

    min_count_to_show_label <- 20

    seurat_cell_types <- plot_cell_types_seurat(cell_type_seurat_counts, "seurat_cluster.projected")
    plots[["seurat_cell_types"]] <- seurat_cell_types
}

output_dir <- file.path(getwd(), "Integration_plots")
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
image_write(images, format = "pdf", path = "RCTD_cell_type_plots.pdf", density = 300)

save(seurat_obj, file = "seurat_obj_integrated.Robj")

# EOF
