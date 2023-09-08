# TOOL spatial-transcriptomics-seurat-spatialdimplot.R: "Seurat v4 -Visualize clusters" (Visualise each cluster separately in SpatialDimPlot.)
# INPUT seurat_spatial_obj_pca.Robj: "Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL spatialdimplot.pdf
# PARAMETER OPTIONAL clusters: "Clusters to plot" TYPE STRING DEFAULT "2, 1" (You can which clusters to plot. If you list multiple clusters, use comma \(,\) as separator.)
# RUNTIME R-4.2.3-single-cell
# TOOLS_BIN ""

# 2022-07-28 IH

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

# Load the R-Seurat-object (called seurat_obj)
load("seurat_spatial_obj_pca.Robj")

# Open the pdf file for plotting
pdf(file = "spatialdimplot.pdf", width = 13, height = 7)

# Run SpatialDimPlot with the chosen clusters
clusters <- trimws(unlist(strsplit(clusters, ",")))
clusters <- strtoi(clusters, base = 0L)
SpatialDimPlot(seurat_obj, cells.highlight = CellsByIdentities(seurat_obj, idents = c(clusters)), facet.highlight = TRUE, ncol = 3)

# close the pdf
dev.off()

## EOF
