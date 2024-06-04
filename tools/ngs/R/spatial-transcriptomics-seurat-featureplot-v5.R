# TOOL spatial-transcriptomics-seurat-featureplot-v5.R: "Seurat v5 -Visualize gene expression" (This tool visualizes gene expression data on top of the tissue histology.)
# INPUT OPTIONAL seurat_object.Robj: "Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL Feature_plot.pdf
# PARAMETER OPTIONAL genes: "Gene name\(s\)" TYPE STRING DEFAULT "Hpca, Ttr" (Name\(s\) of the gene to plot. If you list multiple gene names, use comma \(,\) as separator.)
# PARAMETER OPTIONAL point.size: "Point size in spatial feature plot" TYPE DECIMAL DEFAULT 1.6 (Point size for the plot. Default is 1.6)
# PARAMETER OPTIONAL min_transparency: "Minimum transparency" TYPE DECIMAL DEFAULT 1 (Transparency of the points. Default is 1. Transparency of points with lower expression can be downweighted with lower minimum.)
# PARAMETER OPTIONAL max_transparency: "Maximum transparency" TYPE DECIMAL DEFAULT 1 (Transparency of the points. Default is 1.)
# PARAMETER OPTIONAL color.scale: "Determine color scale based on all genes" TYPE [all:yes, feature:no] DEFAULT feature (Determine whether the color scale is based on all genes or individual genes. By default, the color scale is determined for each gene individually and may differ between genes. If you wish to compare gene expression between different genes, it is useful to set this parameter to "yes" so that the color scale is the same for all genes.)
# RUNTIME R-4.2.3-seurat5
# TOOLS_BIN ""

# 2022-07-25 IH
# 2024-03-21 EP Update to Seurat v5

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(Biobase)

print(package.version("Seurat"))
documentVersion("Seurat", package.version("Seurat"))

# Load seurat object (called seurat_obj)
load("seurat_object.Robj")

# Open the pdf file for plotting
pdf(file = "Feature_plot.pdf", width = 13, height = 7)

# Requested genes
genes <- trimws(unlist(strsplit(genes, ",")))

# Genes in Seurat object
all.genes <- rownames(x = seurat_obj)

# Find case insensitive matches between requested genes and genes in Seurat object
matches <- (unlist(lapply(genes, function(g) grep(paste0("^", g, "$"), all.genes, ignore.case = TRUE, value = TRUE))))

# If none of the requested genes are in the Seurat object
if (identical(matches, character(0))) {
  stop(paste("CHIPSTER-NOTE: ", "None of the genes you requested were not found in the Seurat object."))
}

# Spatial feature plot:
SpatialFeaturePlot(seurat_obj, features = matches, keep.scale = color.scale, pt.size.factor = point.size, alpha = c(min_transparency, max_transparency))
# Close the pdf
dev.off()

## EOF
