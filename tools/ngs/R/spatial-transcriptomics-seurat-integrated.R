# TOOL spatial-transcriptomics-seurat-integrated.R: "Seurat v4 -Integration with scRNA-seq data" (Integrate spatial data with scRNA-seq reference to predict the proportion of different celltypes in the Visium spots.)
# INPUT seurat_obj_subset.Robj: "Seurat object" TYPE GENERIC
# INPUT sc_reference: "Reference scRNA-seq dataset" TYPE GENERIC (Reference single-cell RNA dataset for integration as a Seurat object in .rds format.)
# OUTPUT OPTIONAL seurat_obj_integrated.Robj
# OUTPUT OPTIONAL reference_UMAP_plot.pdf
# RUNTIME R-4.2.3-single-cell
# SLOTS 5
# TOOLS_BIN ""

# 2024-04-24 EP Remove renormalization of subsetted seurat object and add it to the subsetting tool

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

# Load the R-Seurat-object (called seurat_obj)
load("seurat_obj_subset.Robj")

# Load the reference dataset
allen_reference <- readRDS("sc_reference")

# normalise the scRNA-seq reference
allen_reference <- SCTransform(allen_reference, ncells = 3000, verbose = FALSE) %>%
    RunPCA(verbose = FALSE) %>%
    RunUMAP(dims = 1:30)

# Open the pdf file for plotting
pdf(file = "reference_UMAP_plot.pdf", width = 13, height = 7)

# Visualise the reference data
# the annotation is stored in the 'subclass' column of object metadata
DimPlot(allen_reference, group.by = "subclass", label = TRUE)

# Find anchors between a reference and the seurat object
anchors <- FindTransferAnchors(reference = allen_reference, query = seurat_obj, normalization.method = "SCT")

predictions.assay <- TransferData(
    anchorset = anchors, refdata = allen_reference$subclass, prediction.assay = TRUE,
    weight.reduction = seurat_obj[["pca"]], dims = 1:30
)

# add the predictions to the seurat object as a new assay
seurat_obj[["predictions"]] <- predictions.assay

DefaultAssay(seurat_obj) <- "predictions"

# close the pdf
dev.off()

# Save the Robj for the next tool
save(seurat_obj, file = "seurat_obj_integrated.Robj")

# EOF
