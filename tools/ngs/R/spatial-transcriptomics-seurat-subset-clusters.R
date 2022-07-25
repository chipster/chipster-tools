# TOOL spatial-transcriptomics-seurat-subset-clusters.R: "Seurat v4 -Subset out clusters" (Subset out anatomical regions first by taking a subset of clusters.)
# INPUT seurat_spatial_obj_pca.Robj: "Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL seurat_obj_subset1.Robj
# OUTPUT OPTIONAL subset.pdf
# PARAMETER OPTIONAL clusters: "Subset of clusters" TYPE STRING DEFAULT "1,2,3,4,5" (Clusters to subset. If you list multiple clusters, use comma \(,\) as separator, for example "1,2,3,4".)
# PARAMETER OPTIONAL rows: "Choose cell rows" TYPE INTEGER DEFAULT 400 (Check with spatialdimplot which cell rows to remove.)
# PARAMETER OPTIONAL columns: "Choose cell columns" TYPE INTEGER DEFAULT 150 (Check with spatialdimplot which cell columns to remove.)
# RUNTIME R-4.1.0-single-cell

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

# Load the R-Seurat-object (called seurat_obj)
load("seurat_spatial_obj_pca.Robj")

clusters <- trimws(unlist(strsplit(clusters, ",")))
clusters <- strtoi(clusters, base=0L)
#Subset of chosen clusters
seurat_obj_subset <- subset(seurat_obj, idents = c(clusters))

# Open the pdf file for plotting
pdf(file="subset.pdf", width=13, height=7) 

#Visualise which cells to remove with SpatialDimPlot
SpatialDimPlot(seurat_obj_subset,cells.highlight = WhichCells(seurat_obj_subset, expression = slice1_imagerow > rows &  slice1_imagecol > columns))

# close the pdf
dev.off()

# Save the Robj for the next tool
save(seurat_obj, file="seurat_obj_subset1.Robj")

## EOF