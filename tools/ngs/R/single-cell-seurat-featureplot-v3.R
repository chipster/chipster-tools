# TOOL single-cell-seurat-featureplot-v3.R: "Seurat v3 -Visualise features in tSNE plot" (Color single cells on a tSNE dimensional reduction plot according to a feature, i.e. gene expression, PC scores, number of genes detected, etc.)
# INPUT seurat_obj.Robj: "Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL tSNEplot.pdf
# PARAMETER OPTIONAL feature_to_plot: "Feature" TYPE [percent.mt, nCount_RNA, nFeature_RNA, orig.ident, PC1, PC2] DEFAULT percent.mt (Denotes which feature to use for coloring the cells.)
# PARAMETER OPTIONAL point.size: "Point size in tSNE plot" TYPE DECIMAL DEFAULT 1 (Point size for tSNE plot. )
# RUNTIME R-3.6.1


# 09.04.2019 ML
# 09.07.2019 ML Seurat v3

library(Seurat)
library(gplots)

# Load the R-Seurat-object (called seurat_obj)
load("seurat_obj.Robj")

if (exists("data.combined") ){
	seurat_obj <- data.combined
}

# Plot tSNE 
pdf(file="tSNEplot.pdf") 
FeaturePlot(object = seurat_obj, features = feature_to_plot, pt.size = point.size)
dev.off() # close the pdf

# EOF