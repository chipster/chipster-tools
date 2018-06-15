# TOOL single-cell-seurat-pca.R: "Seurat -PCA" (Principal component analysis on the highly variable genes across the single cells.)
# INPUT OPTIONAL seurat_obj.Robj: "Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL PCAplots.pdf
# OUTPUT OPTIONAL PCAgenes.txt
# OUTPUT OPTIONAL seurat_obj_2.Robj
# PARAMETER OPTIONAL num.of.heatmaps: "Number of principal components to plot as heatmaps" TYPE INTEGER DEFAULT 12 (How many principal components to plot as heatmaps.)
# RUNTIME R-3.4.3


# OUTPUT OPTIONAL log.txt

# 2017-06-06 ML
# 2017-07-05 ML split into separate tool
# 2018-01-11 ML update Seurat version to 2.2.0


library(Seurat)
library(dplyr)
library(Matrix)
library(gplots)


# Load the R-Seurat-object (called seurat_obj)
load("seurat_obj.Robj")

# PCA
# The variable genes = genes in seurat_object@var.genes are used as input
seurat_obj <- RunPCA(object = seurat_obj, pc.genes = seurat_obj@var.genes, do.print = TRUE, pcs.print = 5, genes.print = 5)

#PCA genes in txt file
sink("PCAgenes.txt")
PrintPCA(seurat_obj, pcs.print = 1:5, genes.print = 5, use.full = TRUE)
sink()
# PDF
pdf(file="PCAplots.pdf", , width=9, height=12) 
VizPCA(seurat_obj, 1:2)
PCAPlot(seurat_obj, 1, 2)

# Need to check the number of cells at this point.
# cells_left <- seurat_obj@data@Dim[2] 
cells_left <- dim(seurat_obj@raw.data)[2]
if (cells_left > 500) {
	PCHeatmap(seurat_obj, pc.use = 1, cells.use = 100, do.balanced = TRUE)
	PCHeatmap(seurat_obj, pc.use = 1:num.of.heatmaps, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)
}else{
	PCHeatmap(seurat_obj, pc.use = 1, cells.use = cells_left, do.balanced = TRUE)
	PCHeatmap(seurat_obj, pc.use = 1:num.of.heatmaps, cells.use = cells_left, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)
}
# fig.height=12,fig.width=9 
PCElbowPlot(seurat_obj)
dev.off() # close the pdf

# Save the Robj for the next tool
save(seurat_obj, file="seurat_obj_2.Robj")

## EOF
