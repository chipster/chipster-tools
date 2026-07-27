# TOOL spatial-transcriptomics-seurat-integration-with-scRNA-HD-v5.R: "Seurat v5 HD -Integration with scRNA-seq data" (This tool identifies differentially expressed genes between two user defined clusters and visualizes these genes on top of the tissue image.)
# INPUT seurat_obj_clustering.Robj: "Seurat object" TYPE GENERIC
# INPUT seurat_obj_scrna.Robj: "Seurat scRNA data" TYPE GENERIC
# OUTPUT OPTIONAL spatiaaliplotti.pdf
# PARAMETER OPTIONAL label.size: "determine the label size of the plots" TYPE INTEGER DEFAULT 3
# RUNTIME R-4.5.1-seurat5
# SLOTS 6
# TOOLS_BIN ""

resolution = 0.5
dims.reduction = 30
label.size = 2
num_cores = as.numeric(chipster.threads.max)
parallel = TRUE

# 2026-07-15 JV


# Matrix version 1.8.0 needed for zgeMatrix. Trying to install it on the VM Runtime R-4.5.1-seurat5 through github

library(Seurat)
library(ggplot2)
library("Matrix")
library("SeuratWrappers")
library("spacexr")

print("Session info:")
sessionInfo()

load("seurat_obj_clustering.Robj")

spatial_obj <- seurat_obj
rm(seurat_obj)

# For now use sketching, add as a param later
sketch = TRUE
if (sketch) {
    DefaultAssay(spatial_obj) <- "Spatial.008um"

    spatial_obj <- FindVariableFeatures(spatial_obj)

    spatial_obj <- SketchData(
        object = spatial_obj,
        ncells = 50000,
        method = "LeverageScore",
        skeched.assay = "sketch"
    )


    DefaultAssay(spatial_obj) <- "sketch"

    spatial_obj <- ScaleData(spatial_obj)
    spatial_obj <- RunPCA(spatial_obj, assay="sketch", reduction.name = "pca.spatial_obj.sketch", verbose = T)
    spatial_obj <- FindNeighbors(spatial_obj, reduction = "pca.spatial_obj.sketch", dims = 1:50)
    spatial_obj <- RunUMAP(spatial_obj, reduction = "pca.spatial_obj.sketch", reduction.name = "umap.spatial_obj.sketch", return.model = T, dims = 1:50, verbose = T)
}

load("seurat_obj_scrna.Robj")

ref <- seurat_obj
rm(seurat_obj)


Idents(ref) <- "subclass_label"
counts <- ref[["RNA"]]$counts
cluster <- as.factor(ref$subclass_label)
nUMI <- ref$nCount_RNA
levels(cluster) <- gsub("/", "-", levels(cluster))
cluster <- droplevels(cluster)

# create the RCTD reference object
reference <- Reference(counts, cluster, nUMI)

counts_hd <- spatial_obj[["sketch"]]$counts
spatial_obj_cells_hd <- colnames(spatial_obj[["sketch"]])
coords <- GetTissueCoordinates(spatial_obj)[spatial_obj_cells_hd,1:2]

# create the RCTD query object
query <- SpatialRNA(coords, counts_hd, colSums(counts_hd))

# Run RCTD

RCTD <- create.RCTD(query, reference, max_cores = 4)
RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")

spatial_obj <- AddMetaData(spatial_obj, metadata = RCTD@results$results_df)

# Project RCTD labels 

# project RCTD labels from sketched cortical cells to all cortical cells
spatial_obj$first_type <- as.character(spatial_obj$first_type)
spatial_obj$first_type[is.na(spatial_obj$first_type)] <- 'Unknown'
spatial_obj <- ProjectData(
  object = spatial_obj,
  assay = "Spatial.008um",
  full.reduction = "pca.spatial_obj",
  sketched.assay = "sketch",
  sketched.reduction = "pca.spatial_obj.sketch",
  umap.model = "umap.spatial_obj.sketch",
  dims = 1:50,
  refdata = list(full_first_type = "first_type")
)

# Plotting


DefaultAssay(spatial_obj) <- "Spatial.008um"

# we only ran RCTD on the cortical cells
# set labels to all other cells as "Unknown"
spatial_obj[[]][, "full_first_type"] <- "Unknown"
spatial_obj$full_first_type[Cells(spatial_obj)] <- spatial_obj$full_first_type[Cells(spatial_obj)]
Idents(spatial_obj) <- 'full_first_type'

# now we can spatially map the location of any scRNA-seq cell type
# start with Layered (starts with L), excitatory neurons in the cortex
cells <- CellsByIdentities(spatial_obj)
excitatory_names <- sort(grep("^L.* CTX",names(cells),value = TRUE)) # This needs to be customized


pdf(file = "spatiaaliplotti.pdf")
p <- SpatialDimPlot(spatial_obj, cells.highlight = cells[excitatory_names], cols.highlight = c("#FFFF00","grey50"), facet.highlight = T, combine=T, ncol=4)
print(p)

dev.off()

# EOF