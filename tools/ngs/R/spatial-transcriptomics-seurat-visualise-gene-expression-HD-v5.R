# TOOL spatial-transcriptomics-seurat-visualise-gene-expression-HD-v5.R: "Seurat v5 -Visualize gene expression in spatial transcriptomics data" (This tool visualizes the predicted underlying composition of cell types in each spatial spot after integration with scRNA-seq reference data.)
# INPUT seurat_obj_clustering.Robj: "Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL gene_expression_plot.pdf
# PARAMETER OPTIONAL genes: "Genes to plot" TYPE UNCHECKED_STRING DEFAULT "GAPDH, RORB" (If you list multiple genes, please use comma\(s\) \(,\) as a separator, e.g., \"GAPDH\,RORB\".)
# RUNTIME R-4.5.1-seurat5
# SLOTS 4
# TOOLS_BIN ""




# 2026-09-06 JV

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(Biobase)


packageVersion("Seurat")
packageVersion("SeuratObject")
#print(package.version("Seurat"))
#documentVersion("Seurat", package.version("Seurat"))

# Load the R-Seurat-object (called seurat_obj)
load("seurat_obj_clustering.Robj")

genes <- trimws(unlist(strsplit(genes, ",")))

rownames(seurat_obj) <- toupper(rownames(seurat_obj))


colnames(seurat_obj@meta.data)

pdf(file = "gene_expression_plot.pdf", width = 10, height = 10)

p1 <- SpatialFeaturePlot(seurat_obj, features = genes)
p2 <- SpatialDimPlot(seurat_obj, group.by = "seurat_clusters")

print(p1)
print(p2)

dev.off()

# EOF


