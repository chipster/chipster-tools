# TOOL single-cell-seurat-extract-cluster.R: "Seurat v4 BETA -Extract cells in a cluster" (Extract cells in a particular cluster into a new R-object for closer inspection. As input, use R-object after clustering the data. Read the tool manual to see how you can continue the analysis.)
# INPUT seurat_obj.Robj: "Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL log.txt
# OUTPUT seurat_obj_subset.Robj
# OUTPUT OPTIONAL QCplots.pdf 
# PARAMETER cluster.identifier: "Name of the cluster to extract" TYPE STRING DEFAULT "3" (Name of the cluster you wish to extract.)
# RUNTIME R-4.1.0-single-cell


# 03.10.2019 ML
# 2021-10-04 ML Update to Seurat v4

library(Seurat)
library(dplyr)
library(Matrix)
library(gplots)

# Load the R-Seurat-object (called seurat_obj)
load("seurat_obj.Robj")

# QC plots, open pdf:
pdf(file="QCplots.pdf", , width=13, height=7) 

if (exists("data.combined") ){
	seurat_obj <- data.combined
    data.combined <- subset(x = data.combined, idents = cluster.identifier)
    # Save the Robj for the next tool
    save(data.combined, file="seurat_obj_subset.Robj")
}else{
    seurat_obj <- subset(x = seurat_obj, idents = cluster.identifier)
    # Save the Robj for the next tool
    save(seurat_obj, file="seurat_obj_subset.Robj")

    # QC plots, only work for single sample data?
    # % of mito genes (note: they are named either "MT-CO1" or "mt-Co1", have to check both)
    # NOTE: The pattern provided works for human and mouse gene names. You may need to adjust depending on your system of interest
    seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-|^mt-")

    #  Violinplot
    if (sum(is.na(seurat_obj@meta.data$percent.mt)) < 1) {
	    VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
    } else {
    	VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
    }

    # FeatureScatter (v3)
    plot1 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
    plot2 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
    CombinePlots(plots = list(plot1, plot2))
}


# Number of cells:
textplot(paste("\v \v Number of \n \v \v cells: \n \v \v", length(colnames(x=seurat_obj))), halign="center", valign="center", cex=2)

dev.off() # close the pdf

## EOF