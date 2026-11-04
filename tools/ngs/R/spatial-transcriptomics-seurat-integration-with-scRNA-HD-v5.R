# TOOL spatial-transcriptomics-seurat-integration-with-scRNA-HD-v5.R: "Seurat v5 HD -Integration with scRNA-seq data" (This tool identifies differentially expressed genes between two user defined clusters and visualizes these genes on top of the tissue image.)
# INPUT seurat_obj_clustering.Robj: "Seurat object" TYPE GENERIC
# INPUT scRNAseq_ref.Rds: "Seurat scRNA data" TYPE GENERIC
# OUTPUT spatialPlot.pdf
# OUTPUT OPTIONAL seurat_obj_integrated.Robj
# PARAMETER assay: "Assay to use for RCTD" TYPE [Spatial.008um: "Spatial.008um"] DEFAULT Spatial.008um
# PARAMETER known_types: "Cell types to highlight in the plots" TYPE STRING DEFAULT "Astro,Oligo,Endo"
# PARAMETER OPTIONAL celltype_colnames: "Column name of the cell type information in the scRNA-seq Seurat object" TYPE STRING DEFAULT "subclass_label"
# PARAMETER OPTIONAL label.size: "Label size of the plots" TYPE INTEGER DEFAULT 3
# RUNTIME R-4.5.1-visium-hd
# SLOTS 10
# TOOLS_BIN ""

assay <- as.character(assay)
resolution <- 0.5
dims.reduction <- 30
label.size <- 2
num_cores <- as.numeric(chipster.threads.max)
parallel <- TRUE

# 2026-07-15 JV

# Matrix version 1.8.0 needed for zgeMatrix. Trying to install it on the VM Runtime R-4.5.1-seurat5 through github

library(Seurat)
library(ggplot2)
library("dplyr")
library("tidyverse")
library("Matrix")
library("SeuratWrappers")
library("spacexr")


options(future.globals.maxSize = 3000 * 1024^2) # 3000 MiB
#print("Session info:")
#sessionInfo()

load("seurat_obj_clustering.Robj")
# n_cells <- 20000  # target number of cells
# set.seed(42)      # for reproducibility

# cells_to_keep <- sample(colnames(seurat_obj), size = n_cells)
# seurat_obj <- subset(seurat_obj, cells = cells_to_keep)

#seurat_obj

print("Loading ref object")

#load("seurat_obj_scrna.Robj")

# cells_to_keep <- sample(colnames(data), size = n_cells)
# data <- subset(data, cells = cells_to_keep)

ref <- readRDS("scRNAseq_ref.Rds")
#rm(data)

# Basic filter to make sure
# Possibly not needed as it is not in seurat vignette.
#ref <- subset(ref, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)

# # Seurat vignette

# The clustering tool has sketched the 8um assay: the sketch is in the "sketch" assay, its PCA in
# "pca.sketch" and the PCA of all bins in "full.pca.sketch"
# GetTissueCoordinates() below takes the image of the default assay, so set it
DefaultAssay(seurat_obj) <- assay

# Not needed anymore: the clustering tool always sketches the 8um assay, so the input always has
# the "sketch" assay. Kept for reference, e.g. if this tool should accept objects that weren't sketched.
# To use this again:
# - the reduction names must match the ProjectData() call below ("pca.sketch", "full.pca.sketch"),
#   or that call fails
# - FindVariableFeatures() and ScaleData() are needed before SketchData(), otherwise
#   LeverageScore aborts (see the clustering tool)
# - dims should match dims.reduction in ProjectData()
# if (!("sketch" %in% Assays(seurat_obj))) {
#     # Sketch: pick a representative subset of bins (leverage score sampling) into the "sketch" assay

print("Sketching the data")

# Because the data is supposed to be a subset of the whole data

seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- SketchData(
    object = seurat_obj,
    ncells = 50000,
    method = "LeverageScore",
    sketched.assay = "sketch",
    features = VariableFeatures(seurat_obj)
)

# PCA, neighbors and UMAP of the sketched bins only
DefaultAssay(seurat_obj) <- "sketch"
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, assay = "sketch", reduction.name = "pca.sketch", verbose = T)
seurat_obj <- FindNeighbors(seurat_obj, reduction = "pca.sketch", dims = 1:30)
# return.model = TRUE so that ProjectData() could project the UMAP to all bins
seurat_obj <- RunUMAP(
    seurat_obj,
    reduction = "pca.sketch",
    reduction.name = "full.umap.sketch",
    return.model = T,
    dims = 1:30,
    verbose = T
)


gc()

# Here we need to know the exact colname of the celltypes e.g. "celltype"

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

#Idents(ref) <- "subclass_label"

counts <- ref[["RNA"]]$counts
cluster <- as.factor(ref[[celltype_colnames, drop = TRUE]])
nUMI <- ref$nCount_RNA
levels(cluster) <- gsub("/", "-", levels(cluster))
cluster <- droplevels(cluster)

# create the RCTD reference object
reference <- Reference(counts, cluster, nUMI)

counts_hd <- seurat_obj[["sketch"]]$counts
seurat_obj_cells_hd <- colnames(seurat_obj[["sketch"]])
coords_raw <- GetTissueCoordinates(seurat_obj)[seurat_obj_cells_hd, 1:2]


# Coordinates stored in 1 and 2 cols?
#coords <- GetTissueCoordinates(seurat_obj)[seurat_obj_cells_hd,1:2]

# 16um maybe not needed at all, currently wont work due to cell mismatch beteween sketch (8um) and 16um...
#if (assay == "Spatial.016um") {
#    img_16um <- Images(seurat_obj)[[2]]
#    seurat_obj_cells_hd <- colnames(seurat_obj[["Spatial.016um"]])
#    coords_raw <- GetTissueCoordinates(seurat_obj, image = img_16um)[seurat_obj_cells_hd, 1:2]
#} else {
#8um bin
#}

query <- SpatialRNA(coords_raw, counts_hd, colSums(counts_hd))

#SpatialDimPlot(seurat_obj)

print("Running RCTD")
# run RCTD with min of 1 cell
RCTD <- create.RCTD(query, reference, max_cores = num_cores, CELL_MIN_INSTANCE = 25, test_mode = FALSE)
RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
# add results back to Seurat object
seurat_obj <- AddMetaData(seurat_obj, metadata = RCTD@results$results_df)

gc()

# project RCTD labels from sketched cortical cells to all cortical cells
seurat_obj$first_type <- as.character(seurat_obj$first_type)
seurat_obj$first_type[is.na(seurat_obj$first_type)] <- 'Unknown'


# Now this is if we want to re-sketch
# print("Projecting RCTD labels")
# seurat_obj <- ProjectData(
#   object = seurat_obj,
#   assay = "Spatial.008um",
#   full.reduction = "pca.seurat_obj.sketch",
#   sketched.assay = "sketch",
#   sketched.reduction = "pca.seurat_obj.sketch",
#   umap.model = "umap.seurat_obj.sketch",
#   dims = 1:30,
#   refdata = list(full_first_type = "first_type")
# )

# If no resketch, then use this:
seurat_obj <- ProjectData(
    object = seurat_obj,
    assay = "Spatial.008um",
    full.reduction = "full.pca.sketch",
    sketched.assay = "sketch",
    sketched.reduction = "pca.sketch",
    umap.model = "umap.sketch",
    dims = 1:30,
    refdata = list(full_first_type = "first_type")
)

DefaultAssay(seurat_obj) <- assay
Idents(seurat_obj) <- "full_first_type"

# we only ran RCTD on the cortical cells - NOPE ALL
# set labels to all other cells as "Unknown"
#seurat_obj[[]][, "full_first_type"] <- "Unknown"
#seurat_obj$full_first_type[Cells(seurat_obj)] <- seurat_obj$full_first_type[Cells(seurat_obj)]

print("Plotting")

pdf(file = "spatialPlot.pdf", width = 12, height = 8)

p <- SpatialDimPlot(seurat_obj)
print(p)
#, cells.highlight = cells[excitatory_names], cols.highlight = c("#FFFF00","grey50"), facet.highlight = T, combine=T,

# sanity check plot for cortical data (look that the right cells are in right spot physicall)

# seurat_obj[[]][, "full_first_type"] <- "Unknown"
# seurat_obj$full_first_type[Cells(seurat_obj)] <- seurat_obj$full_first_type[Cells(seurat_obj)]
# Idents(seurat_obj) <- "full_first_type"

cells <- CellsByIdentities(seurat_obj)

# excitatory_names <- sort(grep("^L.* CTX",names(cells),value = TRUE))

# p3 <- SpatialDimPlot(seurat_obj, cells.highlight = cells[excitatory_names], cols.highlight = c("#FFFF00","grey50"), facet.highlight = T, combine=T, ncol=4)
# print(p3)

# Claude code

#known_types <- names(sort(table(seurat_obj$full_first_type), decreasing = TRUE))[1:8]
#known_types <- setdiff(known_types, "Unknown")
known_types <- known_types

known_types <- strsplit(known_types, ",")[[1]]
known_types <- trimws(known_types) # Remove leading/trailing whitespace

print(known_types)


p_highlight <- SpatialDimPlot(
    seurat_obj,
    cells.highlight = cells[known_types],
    cols.highlight = c("#FFFF00", "black"),
    facet.highlight = TRUE,
    combine = TRUE,
    ncol = 4,
    image.alpha = 0
)
print(p_highlight)

p_highlight3 <- SpatialDimPlot(
    seurat_obj,
    cells.highlight = cells[known_types],
    cols.highlight = c("#FFFF00", "black"),
    facet.highlight = TRUE,
    combine = TRUE,
    ncol = 4,
    image.alpha = 1
)
print(p_highlight3)

# p_highlight2 <- SpatialDimPlot(
#     seurat_obj,
#     cells.highlight = cells,
#     cols.highlight = c("#FFFF00", "black"),
#     facet.highlight = TRUE,
#     combine = TRUE,
#     ncol = 8,
#     image.alpha = 0.3
# )

p_highlight2 <- SpatialDimPlot(
    seurat_obj,
    cells.highlight = cells[setdiff(names(cells), "NA")],
    cols.highlight = c("#FFFF00", "grey50"),
    facet.highlight = T,
    combine = T,
    ncol = 8
)

print(p_highlight2)

# Seurat vignette code
if ("BANKSY" %in% Assays(seurat_obj)) {
    plot_cell_types <- function(data, label) {
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

    p_banksy <- plot_cell_types(cell_type_banksy_counts, "banksy_cluster")
    print(p_banksy)
}

dev.off()

save(seurat_obj, file = "seurat_obj_integrated.Robj")

# EOF
