# TOOL spatial-transcriptomics-seurat-subset-clusters-v5.R: "Seurat v5 -Subset out anatomical regions based on clusters" (Subset out anatomical regions based on clusters and visualize the result.)
# INPUT seurat_spatial_obj_pca.Robj: "Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL seurat_obj_subset.Robj
# OUTPUT OPTIONAL subset.pdf
# PARAMETER OPTIONAL clusters: "Subset of clusters" TYPE STRING DEFAULT "1,2,3,4,5" (Clusters to subset. If you list multiple clusters, use comma \(,\) as separator, for example "1,2,3,4".)
# PARAMETER OPTIONAL num.features: "Number of variable genes to return in SCTransform" TYPE INTEGER DEFAULT 3000 (After subsetting the data, the subsetted data must be renormalized using SCTransform. Select the number of top variable features, i.e. how many features are returned. For SCTransform, the recommended default is 3000.)
# PARAMETER OPTIONAL PCstocompute: "Number of PCs to compute in PCA" TYPE INTEGER DEFAULT 50 (After renormalizing the data, PCA must be run again. Select the number of PCs to compute in PCA.)
# RUNTIME R-4.2.3-seurat5
# TOOLS_BIN ""

# 2022-08-03 IH
# 2024-03-21 EP Update to Seurat v5

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(Biobase)

print(package.version("Seurat"))
documentVersion("Seurat", package.version("Seurat"))

# Load the R-Seurat-object (called seurat_obj)
load("seurat_spatial_obj_pca.Robj")

clusters <- trimws(unlist(strsplit(clusters, ",")))
clusters <- strtoi(clusters, base = 0L)

# Subset of chosen clusters
seurat_obj <- subset(seurat_obj, idents = c(clusters))

# After subsetting, we renormalize the subsetted spatial data and run PCA again
seurat_obj <- SCTransform(seurat_obj, assay = "Spatial", variable.features.n = num.features, verbose = FALSE)
seurat_obj <- RunPCA(seurat_obj, assay = "SCT", npcs = PCstocompute, verbose = FALSE)

# Open the pdf file for plotting
pdf(file = "subset.pdf", width = 13, height = 7)

SpatialDimPlot(seurat_obj, crop = FALSE, label = TRUE, pt.size.factor = 1, label.size = 3)

# Close the pdf
dev.off()

# Save the Robj for the next tool
save(seurat_obj, file = "seurat_obj_subset.Robj")

## EOF
