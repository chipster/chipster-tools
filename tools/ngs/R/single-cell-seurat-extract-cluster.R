# TOOL single-cell-seurat-extract-cluster.R: "Seurat v3 BETA -Extract cells in a cluster" (Extract cells in a particular cluster into a new R-object for closer inspection.)
# INPUT seurat_obj.Robj: "Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL log.txt
# OUTPUT seurat_obj_subset.Robj
# OUTPUT OPTIONAL QCplots.pdf 
# PARAMETER cluster.identifier: "Name of the cluster to extract" TYPE STRING DEFAULT "3" (Name of the cluster you wish to extract.)
# RUNTIME R-3.6.1


# 03.10.2019 ML

library(Seurat)
library(dplyr)
library(Matrix)
library(gplots)

# Load the R-Seurat-object (called seurat_obj)
load("seurat_obj.Robj")

if (exists("data.combined") ){
	seurat_obj <- data.combined
}

seurat_obj <- subset(x = seurat_obj, idents = cluster.identifier)

# Save the Robj for the next tool
save(seurat_obj, file="seurat_obj_subset.Robj")

# QC
# % of mito genes (note: they are named either "MT-CO1" or "mt-Co1", have to check both)
# NOTE: The pattern provided works for human and mouse gene names. You may need to adjust depending on your system of interest
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-|^mt-")

# pdf plots
pdf(file="QCplots.pdf", , width=13, height=7) 

# Violinplot
if (sum(is.na(seurat_obj@meta.data$percent.mt)) < 1) {
	VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
} else {
	VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
}

# FeatureScatter (v3)
plot1 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

# Number of cells:
textplot(paste("\v \v Number of \n \v \v cells: \n \v \v", length(colnames(x=seurat_obj))), halign="center", valign="center", cex=2)

dev.off() # close the pdf

## EOF