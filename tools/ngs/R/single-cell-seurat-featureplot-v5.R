# TOOL single-cell-seurat-featureplot-v5.R: "Seurat v5 -Visualize features in UMAP plot" (Color single cells on a UMAP dimensional reduction plot according to a feature, i.e. gene expression, PC scores, number of genes detected, etc.)
# INPUT seurat_obj.Robj: "Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL UMAPplot.pdf
# PARAMETER OPTIONAL feature_to_plot: "Feature" TYPE [percent.mt, nCount_RNA, nFeature_RNA] DEFAULT percent.mt (Denotes which feature to use for coloring the cells.)
# PARAMETER OPTIONAL point.size: "Point size in UMAP plot" TYPE DECIMAL DEFAULT 1 (Point size for UMAP plot. )
# PARAMETER OPTIONAL add.labels: "Add labels on top of clusters in plot" TYPE [TRUE: yes, FALSE: no] DEFAULT FALSE (Add cluster number on top of the cluster in UMAP plot.)
# PARAMETER OPTIONAL plotting.order.used: "Plotting order of cells based on expression" TYPE [TRUE:yes, FALSE:no] DEFAULT FALSE (Plot cells in the the order of expression. Can be useful to turn this on if cells expressing given feature are getting buried.)
# RUNTIME R-4.3.2-single-cell
# TOOLS_BIN ""


# 09.04.2019 ML
# 09.07.2019 ML Seurat v3
# 09.09.2019 ML remove optins for orig.ident, PC1, PC2
# 2021-10-04 ML Update to Seurat v4
# 2023-10-25 IH Update to Seurat v5

library(Seurat)
library(gplots)
options(Seurat.object.assay.version = "v5")

# Load the R-Seurat-object (called seurat_obj)
load("seurat_obj.Robj")

if (exists("data.combined")) {
    seurat_obj <- data.combined
}

# Plot tSNE
pdf(file = "UMAPplot.pdf")
FeaturePlot(object = seurat_obj, features = feature_to_plot, pt.size = point.size, label = add.labels, order = as.logical(plotting.order.used))
dev.off() # close the pdf

# EOF
