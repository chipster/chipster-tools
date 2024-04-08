# TOOL spatial-transcriptomics-seurat-featureplot-v5.R: "Seurat v5 -Visualize gene expression" (Visualize molecular data on top of the tissue histology.)
# INPUT OPTIONAL seurat_object.Robj: "Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL Feature_plot.pdf
# PARAMETER OPTIONAL genes: "Gene name\(s\)" TYPE STRING DEFAULT "Hpca, Ttr" (Name\(s\) of the gene to plot. If you list multiple gene names, use comma \(,\) as separator.)
# PARAMETER OPTIONAL point.size: "Point size in spatial feature plot" TYPE DECIMAL DEFAULT 1.6 (Point size for the plot. Default is 1.6)
# PARAMETER OPTIONAL min_transparency: "Minimum transparency" TYPE DECIMAL DEFAULT 1 (Transparency of the points. Default is 1. Transparency of points with lower expression can be downweighted with lower minimum.)
# PARAMETER OPTIONAL max_transparency: "Maximum transparency" TYPE DECIMAL DEFAULT 1 (Transparency of the points. Default is 1.)
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

# Gene expression visualization
genes <- trimws(unlist(strsplit(genes, ",")))

# Sanity check: are all of the requested genes available in the data (one missing allowed)
all.genes <- rownames(x = seurat_obj)
match(genes, all.genes)

# If more than one of the genes is not in the list, print error message:
if (sum(is.na((match(genes, all.genes)))) > 1) {
  not.found <- (genes[is.na(match(genes, all.genes)) == TRUE])
  not.found <- paste(not.found, collapse = ",")
  stop(paste("CHIPSTER-NOTE: ", "The genes you requested were not found in this dataset:", not.found))
}

# Continue even if one gene is missing
if (!all(!is.na(match(genes, all.genes)))) {
  not.found <- genes[is.na(match(genes, all.genes)) == TRUE]
  print(paste("Continuing the visualization without the one gene not found: ", not.found))
  genes <- genes[!is.na(match(genes, all.genes))]
}

# Spatial feature plot:
SpatialFeaturePlot(seurat_obj, features = genes, pt.size.factor = point.size, alpha = c(min_transparency, max_transparency))

# Close the pdf
dev.off()

## EOF
