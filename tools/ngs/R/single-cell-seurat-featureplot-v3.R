# TOOL single-cell-seurat-featureplot-v3.R: "Seurat v3 -Visualise features in UMAP plot" (Color single cells on a UMAP dimensional reduction plot according to a feature, i.e. gene expression, PC scores, number of genes detected, etc.)
# INPUT seurat_obj.Robj: "Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL UMAPplot.pdf
# PARAMETER OPTIONAL feature_to_plot: "Feature" TYPE [percent.mt, nCount_RNA, nFeature_RNA] DEFAULT percent.mt (Denotes which feature to use for coloring the cells.)
# PARAMETER OPTIONAL point.size: "Point size in UMAP plot" TYPE DECIMAL DEFAULT 1 (Point size for UMAP plot. )
# RUNTIME R-3.6.1-single-cell


# 09.04.2019 ML
# 09.07.2019 ML Seurat v3
# 09.09.2019 ML remove optins for orig.ident, PC1, PC2

library(Seurat)
library(gplots)

# Load the R-Seurat-object (called seurat_obj)
load("seurat_obj.Robj")

if (exists("data.combined") ){
	seurat_obj <- data.combined
}

# Plot tSNE 
pdf(file="UMAPplot.pdf") 
FeaturePlot(object = seurat_obj, features = feature_to_plot, pt.size = point.size)
dev.off() # close the pdf

# EOF