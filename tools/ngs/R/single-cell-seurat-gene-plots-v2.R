# TOOL single-cell-seurat-gene-plots-v2.R: "Seurat v2 -Visualize cluster marker genes" (Visualize selected cluster marker genes with violin and feature plot.)
# INPUT seurat_obj.Robj: "Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL log.txt
# OUTPUT OPTIONAL seurat_obj_2.Robj
# OUTPUT OPTIONAL biomarker_plot.pdf
# PARAMETER biomarker: "Gene name" TYPE STRING DEFAULT MS4A1 (Name of the biomarker gene to plot.)
# PARAMETER OPTIONAL point.size: "Point size in tSNE plot" TYPE DECIMAL DEFAULT 1 (Point size for tSNE plot.)
# RUNTIME R-3.4.3


# 13.06.2017 ML
# 2018-01-11 ML Update Seurat version to 2.2.0
# 2019-06-28 EK Add point size parameter for tSNE plot

library(Seurat)
library(dplyr)
library(Matrix)
library(gplots)

# Load the R-Seurat-object (called seurat_obj)
load("seurat_obj.Robj")

# Violin plot:
pdf(file="biomarker_plot.pdf") 
VlnPlot(object = seurat_obj, features.plot = biomarker)
FeaturePlot(object = seurat_obj, features.plot =biomarker, cols.use=c("grey", "blue"), pt.size = point.size, reduction.use = "tsne")
dev.off() # close the pdf

