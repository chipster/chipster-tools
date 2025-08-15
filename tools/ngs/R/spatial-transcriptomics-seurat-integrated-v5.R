# TOOL spatial-transcriptomics-seurat-integrated-v5.R: "Seurat v5 -Integration with scRNA-seq data" (This tool integrates spatial data with scRNA-seq reference data to predict the proportion of different celltypes in the spatial spots.)
# INPUT seurat_obj_subset.Robj: "Spatial Seurat object" TYPE GENERIC
# INPUT sc_reference: "Reference scRNA-seq dataset" TYPE GENERIC (Reference single-cell RNA dataset for integration.)
# OUTPUT OPTIONAL seurat_obj_integrated.Robj
# OUTPUT OPTIONAL reference_UMAP_plot.pdf
# PARAMETER OPTIONAL num.features: "Number of variable genes to return in SCTransform" TYPE INTEGER DEFAULT 3000 (Number of highest variable genes to return in normalization of scRNA-seq reference with SCTransform.)
# PARAMETER OPTIONAL num.cells: "Number of subsampling cells in SCTransform" TYPE INTEGER DEFAULT 5000 (Number of subsampling cells used in learning the noise model in SCTransform for the scRNA-seq reference. Setting this to a smaller number speeds up the computation but might result in loss of performance.)
# PARAMETER OPTIONAL PCstocompute: "Number of PCs to compute" TYPE INTEGER DEFAULT 50 (Number of PCs to compute in performing PCA for scRNA-seq reference.)
# PARAMETER OPTIONAL dims.reduction: "Dimensions of reduction to use" TYPE INTEGER DEFAULT 30 (Dimensions of reduction to use for running UMAP for scRNA-seq reference and performing integration of scRNA-seq data.)
# RUNTIME R-4.2.3-seurat5
# SLOTS 5
# TOOLS_BIN ""

# 2024-03-21 EP Update to Seurat v5

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(Biobase)

print(package.version("Seurat"))
documentVersion("Seurat", package.version("Seurat"))

# Load the R-Seurat-object (called seurat_obj)
load("seurat_obj_subset.Robj")

# Load the reference dataset
allen_reference <- readRDS("sc_reference")

# Normalise the scRNA-seq reference
allen_reference <- SCTransform(allen_reference, ncells = num.cells, variable.features.n = num.features, verbose = FALSE) %>%
    RunPCA(verbose = FALSE) %>%
    RunUMAP(dims = 1:dims.reduction)

# Open the pdf file for plotting
pdf(file = "reference_UMAP_plot.pdf", width = 13, height = 7)

# Visualise the reference data
# The annotation is stored in the 'subclass' column of object metadata
DimPlot(allen_reference, group.by = "subclass", label = TRUE)

# Find anchors between a reference and the seurat object
anchors <- FindTransferAnchors(reference = allen_reference, query = seurat_obj, normalization.method = "SCT")

predictions.assay <- TransferData(
    anchorset = anchors, refdata = allen_reference$subclass, prediction.assay = TRUE,
    weight.reduction = seurat_obj[["pca"]], dims = 1:dims.reduction
)

# Add the predictions to the seurat object as a new assay
seurat_obj[["predictions"]] <- predictions.assay

DefaultAssay(seurat_obj) <- "predictions"

# Close the pdf
dev.off()

# Save the Robj for the next tool
save(seurat_obj, file = "seurat_obj_integrated.Robj")

# EOF
