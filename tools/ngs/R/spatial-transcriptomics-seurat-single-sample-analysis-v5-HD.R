# TOOL spatial-transcriptomics-seurat-single-sample-analysis-v5-HD.R: "Seurat v5 HD -Normalization and PCA" (This tool performs normalization and PCA for a single Seurat object. It first normalizes data with SCTransform and detects highly variable genes. Then, it performs principal component analysis on the highly variable genes detected by SCTransform.)
# INPUT seurat_object.Robj: "Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL seurat_spatial_obj_pca.Robj
# OUTPUT OPTIONAL PCAloadings.txt
# OUTPUT OPTIONAL PCAplots.pdf
# PARAMETER OPTIONAL num.features: "Number of variable genes to return in SCTransform" TYPE INTEGER DEFAULT 3000 (Number of highest variable genes to return in SCTransform.)
# PARAMETER OPTIONAL PCstocompute: "Number of PCs to compute" TYPE INTEGER DEFAULT 50 (Number of PCs to compute in PCA.)
# PARAMETER OPTIONAL loadings: "Print loadings in a file" TYPE [TRUE: yes, FALSE: no] DEFAULT FALSE (Print the PC loadings in a txt file.)
# PARAMETER OPTIONAL num.of.genes.loadings: "Number of genes to list in the loadings file" TYPE INTEGER DEFAULT 5 (How many genes to list in the loadings txt file.)
# RUNTIME R-4.2.3-seurat5
# SLOTS 3
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

# Bins / Assay
# PARAMETER OPTIONAL bin_sizes: "Bin sizes" TYPE STRING DEFAULT "8, 16" (List here the bin sizes, separated by comma)
# bin_sizes <- as.numeric(trimws(unlist(strsplit(bin_sizes, ","))))
# print(bin_sizes)

# NEEDS THE LOOP STILL.
assay_names <- Assays(seurat_obj) # [1] "Spatial.008um" "Spatial.016um"
DefaultAssay(seurat_obj) <- assay_names[1] # "Spatial.008um"
assay_bin <- assay_names[1]
nCount_bin <- paste("nCount_",assay_names[1], sep="")
nFeature_bin <- paste("nFeature_",assay_names[1], sep="")
just.bin <- sub("Spatial\\.", "", assay_names[1])
pca_bin <- paste("pca.", just.bin, sep="")

# SCTransform normalizes the data, detects high-variance features, and stores the data in the SCT assay.
# "--we generally don't recommend using sctransform on spatial data, instead you can run through the older version of normalization, scaling, and find variable features."
# https://github.com/satijalab/seurat/issues/9456
# seurat_obj <- SCTransform(seurat_obj, assay = assay_bin, variable.features.n = num.features, verbose = FALSE)
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, reduction.name = pca_bin, npcs = PCstocompute, verbose = FALSE )
# seurat_obj <- RunPCA(seurat_obj, assay = "SCT", npcs = PCstocompute, verbose = FALSE)

# The following makes sure that all "umi.assay" slots are called "Spatial" after SCTransform
# This is not an issue with one sample. It is an issue, when using this tool after subsetting out anatomical 
# regions based on clusters from a Seurat object that has been merged/integrated. Issue is explained by StepahieHowe 
# in issue https://github.com/satijalab/seurat/issues/8235. If the user is trying to normalize a Seurat object 
# that has been merged or integrated, there will be an error in the "Seurat v5 -Integration with scRNA-seq
# data" tool. This solution does not solve issues with FindSpatiallyVariableFeatures() but it allows to
# perform integration with scRNA-seq data without any errors. FindSpatiallyVariableFeatures() must still be
# run for each object individually in the next tool because otherwise it will result in 'please provide the 
# same number of observations as spatial locations' error as explained by tingchiafelix
# print(length(seurat_obj@assays$SCT@SCTModel.list))
# for (i in 1:length(seurat_obj@assays$SCT@SCTModel.list)) {
#     slot(object = seurat_obj@assays$SCT@SCTModel.list[[i]], name="umi.assay")<-"Spatial"
# }
    
# You can use this to check that all "umi.assay" slots are called "Spatial"
print(SCTResults(object=seurat_obj, slot="umi.assay"))


# PCA genes in txt file
if (loadings == TRUE) {
    sink("PCAloadings.txt")
    print(seurat_obj[["pca"]], dims = 1:PCstocompute, nfeatures = num.of.genes.loadings)
    sink()
}

# PDF plots
pdf(file = "PCAplots.pdf", , width = 9, height = 12)
ElbowPlot(seurat_obj, ndims = PCstocompute) + ggtitle("Amount of variation in the data explained by each PC")
VizDimLoadings(seurat_obj, dims = 1:1, reduction = pca_bin) + ggtitle("Top 30 genes associated with PC 1")
VizDimLoadings(seurat_obj, dims = 2:2, reduction = pca_bin) + ggtitle("Top 30 genes associated with PC 2")
dev.off() # close the pdf

# Save the Robj for the next tool
save(seurat_obj, file = "seurat_spatial_obj_pca.Robj")

# EOF
