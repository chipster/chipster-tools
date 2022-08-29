# TOOL spatial-transcriptomics-seurat-subset-clusters.R: "Seurat v4 -Subset out anatomical regions based on clusters" (Subset out anatomical regions based on clusters and visualise the result.)
# INPUT seurat_spatial_obj_pca.Robj: "Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL seurat_obj_subset.Robj
# OUTPUT OPTIONAL subset.pdf
# PARAMETER OPTIONAL clusters: "Subset of clusters" TYPE STRING DEFAULT "1,2,3,4,5" (Clusters to subset. If you list multiple clusters, use comma \(,\) as separator, for example "1,2,3,4".)
# RUNTIME R-4.1.0-single-cell

# 2022-08-03 IH

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

# Load the R-Seurat-object (called seurat_obj)
load("seurat_spatial_obj_pca.Robj")

clusters <- trimws(unlist(strsplit(clusters, ",")))
clusters <- strtoi(clusters, base=0L)

#Subset of chosen clusters
seurat_obj <- subset(seurat_obj, idents = c(clusters))

# Open the pdf file for plotting
pdf(file="subset.pdf", width=13, height=7) 

SpatialDimPlot(seurat_obj, crop = FALSE, label = TRUE, pt.size.factor = 1, label.size = 3)

# close the pdf
dev.off()

# Save the Robj for the next tool
save(seurat_obj, file="seurat_obj_subset.Robj")

## EOF