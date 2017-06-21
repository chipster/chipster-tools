# TOOL single-cell-seurat-gene-plots.R: "BETA Seurat -Visualize biomarkers" (Visualize selected biomarkers with violin and feature plot.)
# INPUT seurat_obj.Robj: "Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL log.txt
# OUTPUT OPTIONAL seurat_obj_2.Robj
# OUTPUT OPTIONAL biomarker_plot.pdf
# PARAMETER biomarker: "Gene name" TYPE STRING DEFAULT MS4A1 (Name of the biomarker gene to plot.)
# RUNTIME R-3.3.2


# 13.06.2017 ML

library(Seurat)
library(dplyr)
library(Matrix)
library(gplots)

# Load the R-Seurat-object (called seurat_obj)
load("seurat_obj.Robj")

# Violin plot:
pdf(file="biomarker_plot.pdf") 
VlnPlot(seurat_obj, biomarker)
FeaturePlot(seurat_obj, biomarker, cols.use=c("grey", "blue"), pt.size=5)
dev.off() # close the pdf

