# TOOL single-cell-seurat-CCA.R: "Seurat v4 -Combine two samples" (This tool can be used to integrate data and combine two Seurat objects for later joined analysis. The two objects need to be named when created in the Seurat setup tool.) 
# INPUT OPTIONAL seurat_obj1.Robj: "Seurat object 1" TYPE GENERIC
# INPUT OPTIONAL seurat_obj2.Robj: "Seurat object 2" TYPE GENERIC
# OUTPUT OPTIONAL CCAplot.pdf
# OUTPUT seurat_obj_combined.Robj
# PARAMETER OPTIONAL CCstocompute: "Number of CCs to use in the neighbor search" TYPE INTEGER DEFAULT 20 (Which dimensions to use from the CCA to specify the neighbor search space. The neighbors are used to determine the anchors for the alignment.)
# PARAMETER OPTIONAL PCstocompute: "Number of PCs to use in the anchor weighting" TYPE INTEGER DEFAULT 20 (Number of PCs to use in the anchor weighting procedure. The anchors and their weights are used to compute the correction vectors, which allow the datasets to be integrated.)
# IMAGE comp-20.04-r-base
# RUNTIME R-4.1.0-single-cell
# SLOTS 3


# SLOTS = 3: when testing at VM this tool required 18.8G)
# Not used atm: PARAMETER OPTIONAL CCstovisualise: "How many CCs to visualise as heatmaps" TYPE INTEGER DEFAULT 9 (How many canonical components you want to visualise as heatmaps.)
# Not used atm: PARAMETER OPTIONAL num.features: "Number of variable features to return" TYPE INTEGER DEFAULT 2000 (How many features returned per dataset.)


# 2018-08-05 ML
# 2018-10-03 ML Add sample identifiers to cell barcodes (fix problem with same cell barcodes in two samples)
# 09.07.2019 ML Seurat v3
# 2021-10-04 ML Update to Seurat v4

library(Seurat)

# Load the R-Seurat-objects (called seurat_obj -that's why we need to rename them here)
load("seurat_obj1.Robj")
seurat_obj1 <- seurat_obj
load("seurat_obj2.Robj")
seurat_obj2 <- seurat_obj

# Perform integration: 
# 1. identify anchors using the FindIntegrationAnchors function
data.anchors <- FindIntegrationAnchors(object.list = list(seurat_obj1, seurat_obj2), dims = 1:CCstocompute) # dims = Which dimensions to use from the CCA to specify the neighbor search space
# 2. use these anchors to integrate the two datasets together with IntegrateData.
data.combined <- IntegrateData(anchorset = data.anchors, dims = 1:PCstocompute) # dims = Number of PCs to use in the weighting procedure

DefaultAssay(data.combined) <- "integrated"

# Note: these steps are now done twice?
data.combined <- ScaleData(data.combined, verbose = FALSE)  
data.combined <- RunPCA(data.combined, npcs = 30, verbose = FALSE)

# No plots, makes it confusing (no CCA option, only PCA, tSNE and UMAP)
# Plots:
# pdf(file="CCAplot.pdf", , width=13, height=7)  # open pdf
# DimPlot(data.combined, reduction = "pca")

# How to decide number of PCAs/CCAs??
# ElbowPlot(data.combined)

# dev.off() # close the pdf

# Save the Robj for the next tool
save(data.combined, file="seurat_obj_combined.Robj")

## EOF



