# TOOL spatial-transcriptomics-seurat-subset-clusters-v5.R: "Seurat v5 -Subset out anatomical regions based on clusters" (This tool subsets out anatomical regions based on clusters and visualizes the subsetted clusters on top of the tissue histology.)
# INPUT seurat_spatial_obj_pca.Robj: "Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL seurat_obj_subset.Robj
# OUTPUT OPTIONAL subset.pdf
# PARAMETER OPTIONAL clusters: "Subset of clusters" TYPE STRING DEFAULT "1,2,3,4,5" (Clusters to subset. If you list multiple clusters, use comma \(,\) as separator, for example "1,2,3,4".)
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

# Open the pdf file for plotting
pdf(file = "subset.pdf", width = 13, height = 7)

SpatialDimPlot(seurat_obj, crop = FALSE, label = TRUE, pt.size.factor = 1, label.size = 3)

# Close the pdf
dev.off()

# Save the Robj for the next tool
save(seurat_obj, file = "seurat_obj_subset.Robj")

## EOF
