# TOOL single-cell-seurat-gene-plots.R: "Seurat v4 -Visualize genes" (Visualize for example selected cluster marker genes with violin and feature plot.)
# INPUT seurat_obj.Robj: "Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL log.txt
# OUTPUT OPTIONAL seurat_obj_2.Robj
# OUTPUT OPTIONAL biomarker_plot.pdf
# PARAMETER biomarker: "Gene name\(s\)" TYPE STRING DEFAULT "MS4A1, LYZ" (Name\(s\) of the biomarker gene to plot. If you list multiple gene names, use comma \(,\) as separator.)
# PARAMETER OPTIONAL point.size: "Point size in cluster plot" TYPE DECIMAL DEFAULT 1 (Point size for tSNE and UMAP plots.)
# PARAMETER OPTIONAL reduction.method: "Visualisation with tSNE, UMAP or PCA" TYPE [umap:UMAP, tsne:tSNE, pca:PCA] DEFAULT umap (Which dimensionality reduction plot to use.)
# IMAGE comp-20.04-r-base
# RUNTIME R-4.1.0-single-cell



# 13.06.2017 ML
# 2018-01-11 ML update Seurat version to 2.2.0
# 2019-06-28 EK Add point size parameter for tSNE plot
# 2019-06-13 ML Seurat v3
# 2020-06-18 ML Add ridge plot
# 2021-10-04 ML Update to Seurat v4

# for UMAP:
library(reticulate)
use_python("/opt/chipster/tools/miniconda3/envs/chipster_tools/bin/python")


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

# Sanity check: are the requested genes available in the data:
all.genes <- rownames(x = seurat_obj)
match(biomarker, all.genes)
# if one of the genes is not in the list, print error message:
if (!all(!is.na(match(biomarker, all.genes)))) { 
  not.found <- biomarker[is.na(match(biomarker, all.genes))==TRUE]
 #  print(paste("The gene you requested was not found in this dataset:", not.found))
  stop(paste('CHIPSTER-NOTE: ', "The gene you requested was not found in this dataset:", not.found))
  }

# open pdf
pdf(file="biomarker_plot.pdf", width=12, height=12) 

# Violin plot:
VlnPlot(seurat_obj, features = biomarker)

# Feature plot:
FeaturePlot(seurat_obj, features = biomarker, pt.size=point.size, reduction=reduction.method) 	

# Ridge plot:
RidgePlot(seurat_obj, features = biomarker, ncol = 2)

# close the pdf
dev.off() 

# EOF