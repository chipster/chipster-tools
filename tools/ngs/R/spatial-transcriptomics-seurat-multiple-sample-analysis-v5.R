# TOOL spatial-transcriptomics-seurat-multiple-sample-analysis-v5.R: "Seurat v5 -Combine multiple samples, normalization, PCA, clustering and visualization" (This tool performs a combined analysis by merging or integrating multiple data sets across sections. The combined analysis includes similar steps as the single sample analysis which includes normalization, PCA, clustering and visualization. Use this tool instead of the single sample analysis tool when you want to analyze multiple samples at the same time.)
# INPUT samples{...}.Robj: "Samples to combine" TYPE GENERIC
# OUTPUT OPTIONAL seurat_obj_multiple.Robj
# OUTPUT OPTIONAL merged_plot.pdf
# OUTPUT OPTIONAL integrated_plot.pdf
# PARAMETER OPTIONAL num.features: "Number of variable genes to return in SCTransform" TYPE INTEGER DEFAULT 3000 (Number of features to select as top variable features, i.e. how many features returned. For SCTransform, the recommended default is 3000.)
# PARAMETER OPTIONAL method: "Combining method" TYPE [merge: Merge, integration: Integration] DEFAULT merge (User can choose to merge or integrate the samples.)
# PARAMETER OPTIONAL PCstocompute: "Number of PCs to compute" TYPE INTEGER DEFAULT 50 (Number of PCs to compute in PCA.)
# PARAMETER OPTIONAL loadings: "Print loadings in a file" TYPE [TRUE: yes, FALSE: no] DEFAULT FALSE (Print the PC loadings to a txt file.)
# PARAMETER OPTIONAL num.of.genes.loadings: "Number of genes to list in the loadings file" TYPE INTEGER DEFAULT 5 (How many genes to list in the loadings txt file.)
# PARAMETER OPTIONAL dims.reduction: "Dimensions of reduction to use as input" TYPE INTEGER DEFAULT 30 (Number of dimensions of reduction to use for clustering and UMAP. If integration is used, this is also the number of dimensions of reduction to use for integration.)
# PARAMETER OPTIONAL res: "Resolution for granularity for clustering" TYPE DECIMAL DEFAULT 0.8 (Resolution parameter that sets the granularity of the clustering. Increased values lead to greater number of clusters.)
# RUNTIME R-4.2.3-seurat5
# SLOTS 3
# TOOLS_BIN ""

# 2022-07-21 IH
# 2022-10-13 ML nfeatures = 3000 fix
# 2024-03-21 EP Update to Seurat v5 (and add SCTransform to this tool)

library(Seurat)
library(gplots)
library(ggplot2)
library(patchwork)
library(dplyr)
library(Biobase)
require(cowplot)
library(Matrix)

print(package.version("Seurat"))
documentVersion("Seurat", package.version("Seurat"))

# Load the R-Seurat-objects (called seurat_obj -that's why we need to rename them here)
input <- read.table("chipster-inputs.tsv", header = FALSE, sep = "\t")
for (i in 1:nrow(input)) {
    load(input[i, 1])
    name <- paste("seurat_obj_", i, sep = "")
    assign(name, seurat_obj)
}
seurat_objects <- as.list(mget(objects(pattern = "seurat_obj_")))

# Check the number of samples
if (length(seurat_objects) < 2) {
    stop(paste("CHIPSTER-NOTE: ", "It seems you don't have multiple samples. Please check your input files."))
}

# Merge multiple slices
if (length(seurat_objects) == 2) {
    objects_combined <- merge(seurat_objects[[1]], seurat_objects[[2]])
} else {
    objects_combined <- merge(seurat_objects[[1]], y = c(seurat_objects[c(2, (length(seurat_objects)))]))
}

# Rename combined Robj
seurat_obj <- objects_combined
print(seurat_obj)

# Split datasets and process first without integration
#seurat_obj[["Spatial"]] <- split(seurat_obj[["Spatial"]], f = seurat_obj$orig.ident)

# SCTransform normalizes the data, detects high-variance features, and stores the data in the SCT assay
seurat_obj <- SCTransform(seurat_obj, assay = "Spatial", variable.features.n = num.features, verbose = FALSE)

seurat_obj <- RunPCA(seurat_obj, assay = "SCT", npcs = PCstocompute, verbose = FALSE)

# PCA genes in txt file
if (loadings == TRUE) {
    sink("PCAloadings.txt")
    print(seurat_obj[["pca"]], dims = 1:PCstocompute, nfeatures = num.of.genes.loadings)
    sink()
}

if (method == "merge") {
    seurat_obj <- FindNeighbors(seurat_obj, dims = 1:dims.reduction, verbose=FALSE)
    seurat_obj <- FindClusters(seurat_obj, resolution = res, verbose = FALSE)
    seurat_obj <- RunUMAP(seurat_obj, dims = 1:dims.reduction, verbose = FALSE)
    
    # Visualization
    pdf(file = "merged_plot.pdf", , width = 9, height = 12)
    print(DimPlot(seurat_obj, reduction = "umap", group.by = "ident"))
    print(DimPlot(seurat_obj, reduction = "umap", group.by = "orig.ident"))
    print(SpatialDimPlot(seurat_obj, label = TRUE, label.size = 3))

    # Close pdf
    dev.off()
}

if (method == "integration") {

    # Name of new reduction
    new.reduction = "integrated.cca"

    # Integration
    seurat_obj <- IntegrateLayers(object = seurat_obj, method = "CCAIntegration", normalization.method = "SCT",
    orig.reduction = "pca", new.reduction = new.reduction, dims = 1:dims.reduction , assay = "SCT", verbose = FALSE)

    seurat_obj <- FindNeighbors(seurat_obj, reduction = new.reduction, dims = 1:dims.reduction, verbose = FALSE)
    seurat_obj <- FindClusters(seurat_obj, resolution = res, verbose = FALSE)
    seurat_obj <- RunUMAP(seurat_obj, reduction = new.reduction, dims = 1:dims.reduction, verbose = FALSE)

    # Visualization
    pdf(file = "integrated_plot.pdf", width = 9, height = 12) 

    print(DimPlot(seurat_obj, reduction = "umap", group.by = "ident"))
    print(DimPlot(seurat_obj, reduction = "umap", group.by = "orig.ident"))
    print(SpatialDimPlot(seurat_obj, label = TRUE, label.size = 3))

    # Close the pdf
    dev.off()
   
}

# Save the combined Robj for the next tool
save(seurat_obj, file = "seurat_obj_multiple.Robj")

# EOF
