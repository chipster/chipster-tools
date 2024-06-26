# TOOL spatial-transcriptomics-seurat-top-expressed-genes-v5.R: "Seurat v5 -Identify top expressed genes" (This tool identifies top expressed genes in the data and visualizes them in a box plot.)
# INPUT OPTIONAL seurat_object.Robj: "Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL top_genes.pdf
# RUNTIME R-4.2.3-seurat5
# TOOLS_BIN ""

# 2022-07-22 IH
# 2024-03-21 EP Update to Seurat v5

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(Matrix)
library(Biobase)

print(package.version("Seurat"))
documentVersion("Seurat", package.version("Seurat"))

# Load the R-Seurat-object (called seurat_obj)
load("seurat_object.Robj")

# Open the pdf file for plotting
pdf(file = "top_genes.pdf", width = 13, height = 7)

# Identify top expressed genes
# code from  https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/seurat/seurat_07_spatial.html
C <- LayerData(seurat_obj, assay = "Spatial", layer = "counts")
C@x <- C@x / rep.int(colSums(C), diff(C@p))
most_expressed <- order(Matrix::rowSums(C), decreasing = T)[20:1]
boxplot(as.matrix(t(C[most_expressed, ])),
    cex = 0.1, las = 1, xlab = "% total count per cell",
    col = (scales::hue_pal())(20)[20:1], horizontal = TRUE
)

# close the pdf
dev.off()

# EOF
