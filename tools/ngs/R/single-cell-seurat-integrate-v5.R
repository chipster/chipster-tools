# TOOL single-cell-seurat-integrate-v5.R: "Seurat v5 -Integrate multiple samples" (This tool integrates multiple samples to match shared cell types and states across dataset. The samples \/R-objects need to be named when created in the Seurat Setup tool.)
# INPUT OPTIONAL seurat_obj_combined.Robj: "Combined seurat object to integrate" TYPE GENERIC
# OUTPUT seurat_obj_integrated.Robj
# PARAMETER OPTIONAL normalisation.method: "Normalisation method used previously" TYPE [LogNormalize:"Global scaling normalization", SCT:"SCTransform"] DEFAULT LogNormalize (Which normalisation method was used in preprocessing, Global scaling normalization \(default, NormalizeData function used\) or SCTransform.)
# PARAMETER OPTIONAL anchor.identification.method: "Anchor identification method" TYPE [cca:CCA, rpca:RPCA] DEFAULT cca (Which anchor identification method to use. By default, canonical correlation analysis CCA is used, but user can also decide to use the faster and more conservative reciprocal PCA approach. Check from the manual in which cases this option is recommended.)
# PARAMETER OPTIONAL CCstocompute: "Number of CCs to use in the neighbor search" TYPE INTEGER DEFAULT 30 (Which dimensions to use from the CCA to specify the neighbor search space. The neighbors are used to determine the anchors for the alignment.)
# PARAMETER OPTIONAL PCstocompute: "Number of PCs to use in the anchor weighting" TYPE INTEGER DEFAULT 30 (Number of PCs to use in the anchor weighting procedure. The anchors and their weights are used to compute the correction vectors, which allow the datasets to be integrated.)
# PARAMETER OPTIONAL ref.sample.names: "Samples to use as references" TYPE STRING DEFAULT "No references selected" (Names of the sample or samples you wish to use as references in integration, separated by comma. If you are integrating several large datasets, the tool might run out of memory. Choosing to use only some of them as references makes the integration more memory efficient and faster. Please note that the sample names here are case sensitive, so check how you typed the names of the samples when running the setup tool.)
# RUNTIME R-4.3.1-single-cell
# SLOTS 2
# TOOLS_BIN ""


# 2021-12-30 ML
# 2022-02-17 EK increased slots to 4
# 2022-04-19 ML increased slots to 5
# 2022-05-04 ML add RPCA option for anchor identification
# 2022-05-05 ML Rewrite the code, add option to use only part of samples as references
# 2023-02-03 ML Add 5 slots
# 2023-04-06 LG Remove 5 slots
# 2023-02-01 ML Return to the original 3 slots
# 2023-11-29 IH Update to Seurat v5


library(BPCells)
library(Seurat)
options(Seurat.object.assay.version = "v5")

# Load the R-Seurat-object
load("seurat_obj_combined.Robj")


#seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
#seurat_obj <- RunPCA(seurat_obj, npcs = 30, verbose = FALSE)
# Compute UMAP without integration
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30, reduction = "pca")
seurat_obj <- FindClusters(seurat_obj, resolution = 2, cluster.name = "unintegrated_clusters")
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")


data.combined <- IntegrateLayers(object = seurat_obj, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca",
    verbose = FALSE)
data.combined[["RNA"]] <- JoinLayers(data.combined[["RNA"]])

data.combined <- FindNeighbors(data.combined, reduction = "integrated.cca", dims = 1:30)
data.combined <- FindClusters(data.combined, resolution = 1)

data.combined <- RunUMAP(data.combined, dims = 1:30, reduction = "integrated.cca")


# Save the Robj for the next tool
save(data.combined, file = "seurat_obj_integrated.Robj")

## EOF
