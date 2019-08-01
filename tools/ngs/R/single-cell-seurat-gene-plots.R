# TOOL single-cell-seurat-gene-plots.R: "Seurat -Visualize cluster marker genes" (Visualize selected cluster marker genes with violin and feature plot.)
# INPUT seurat_obj.Robj: "Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL log.txt
# OUTPUT OPTIONAL seurat_obj_2.Robj
# OUTPUT OPTIONAL biomarker_plot.pdf
# PARAMETER biomarker: "Gene name\(s\)" TYPE STRING DEFAULT MS4A1 (Name\(s\) of the biomarker gene to plot. If you list multiple gene names, use comma \(,\) as separator.)
# PARAMETER OPTIONAL point.size: "Point size in tSNE plot" TYPE DECIMAL DEFAULT 1 (Point size for tSNE plot.)
# RUNTIME R-3.4.3


# 13.06.2017 ML
# 2018-01-11 ML update Seurat version to 2.2.0
# 2019-06-28 EK Add point size parameter for tSNE plot
# 2019-06-13 ML Seurat v3

library(Seurat)
library(dplyr)
library(Matrix)
library(gplots)

# Load the R-Seurat-object (called seurat_obj)
load("seurat_obj.Robj")

if (exists("data.combined") ){
	seurat_obj <- data.combined
}

# If multiple genes are listed: (separate words from "," and remove whitespace)
if(length(grep(",", biomarker)) != 0) {
   biomarker <- trimws(unlist(strsplit(biomarker, ",")))
}

# Violin plot:
pdf(file="biomarker_plot.pdf") 
VlnPlot(seurat_obj, features = biomarker)
FeaturePlot(seurat_obj, features = biomarker, pt.size=point.size) 	
dev.off() # close the pdf

# EOF