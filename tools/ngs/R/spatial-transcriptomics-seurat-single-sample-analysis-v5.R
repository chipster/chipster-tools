# TOOL spatial-transcriptomics-seurat-single-sample-analysis-v5.R: "Seurat v5 -Normalization and PCA" (This tool performs normalization and PCA for a single Seurat object. It first normalizes data with SCTransform and detects highly variable genes. Then, it performs principal component analysis on the highly variable genes detected by SCTransform.)
# INPUT OPTIONAL seurat_object.Robj: "Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL seurat_spatial_obj_pca.Robj
# OUTPUT OPTIONAL PCAloadings.txt
# OUTPUT OPTIONAL PCAplots.pdf
# PARAMETER OPTIONAL num.features: "Number of variable genes to return in SCTransform" TYPE INTEGER DEFAULT 3000 (Number of highest variable genes to return in SCTransform.)
# PARAMETER OPTIONAL PCstocompute: "Number of PCs to compute" TYPE INTEGER DEFAULT 50 (Number of PCs to compute in PCA.)
# PARAMETER OPTIONAL loadings: "Print loadings in a file" TYPE [TRUE: yes, FALSE: no] DEFAULT FALSE (Print the PC loadings in a txt file.)
# PARAMETER OPTIONAL num.of.genes.loadings: "Number of genes to list in the loadings file" TYPE INTEGER DEFAULT 5 (How many genes to list in the loadings txt file.)
# RUNTIME R-4.2.3-seurat5
# TOOLS_BIN ""

# 2022-07-27 IH
# 2022-10-20 ML Add UMAP plot with sample names
# 2024-03-21 EP Update to Seurat v5 (and add SCTransform to this tool)

library(Seurat)
library(dplyr)
library(Matrix)
library(gplots)
library(ggplot2)
library(patchwork)
library(Biobase)

print(package.version("Seurat"))
documentVersion("Seurat", package.version("Seurat"))

# Load seurat object (called seurat_obj)
load("seurat_object.Robj")

# SCTransform normalizes the data, detects high-variance features, and stores the data in the SCT assay.
seurat_obj <- SCTransform(seurat_obj, assay = "Spatial", variable.features.n = num.features, verbose = FALSE)

# The following makes sure that all "umi.assay" slots are called "Spatial" after SCTransform
# This is not an issue with one sample. It is an issue, when using this tool after subsetting out anatomical 
# regions based on clusters from a Seurat object that has been merged/integrated. Issue is explained by StepahieHowe 
# in issue https://github.com/satijalab/seurat/issues/8235. If the user is trying to normalize a Seurat object 
# that has been merged or integrated, there will be an error in the "Seurat v5 -Integration with scRNA-seq
# data" tool. This solution does not solve issues with FindSpatiallyVariableFeatures() but it allows to
# perform integration with scRNA-seq data without any errors. FindSpatiallyVariableFeatures() must still be
# run for each object individually in the next tool because otherwise it will result in 'please provide the 
# same number of observations as spatial locations' error as explained by tingchiafelix
print(length(seurat_obj@assays$SCT@SCTModel.list))
for (i in 1:length(seurat_obj@assays$SCT@SCTModel.list)) {
    slot(object = seurat_obj@assays$SCT@SCTModel.list[[i]], name="umi.assay")<-"Spatial"
}
    
# You can use this to check that all "umi.assay" slots are called "Spatial"
print(SCTResults(object=seurat_obj, slot="umi.assay"))

seurat_obj <- RunPCA(seurat_obj, assay = "SCT", npcs = PCstocompute, verbose = FALSE)

# PCA genes in txt file
if (loadings == TRUE) {
    sink("PCAloadings.txt")
    print(seurat_obj[["pca"]], dims = 1:PCstocompute, nfeatures = num.of.genes.loadings)
    sink()
}

# PDF plots
pdf(file = "PCAplots.pdf", , width = 9, height = 12)
ElbowPlot(seurat_obj, ndims = PCstocompute) + ggtitle("Amount of variation in the data explained by each PC")
VizDimLoadings(seurat_obj, dims = 1:1, reduction = "pca") + ggtitle("Top 30 genes associated with PC 1")
VizDimLoadings(seurat_obj, dims = 2:2, reduction = "pca") + ggtitle("Top 30 genes associated with PC 2")
dev.off() # close the pdf

# Save the Robj for the next tool
save(seurat_obj, file = "seurat_spatial_obj_pca.Robj")

# EOF
