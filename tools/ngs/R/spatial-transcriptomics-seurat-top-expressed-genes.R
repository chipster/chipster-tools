# TOOL spatial-transcriptomics-seurat-top-expressed-genes.R: "Seurat v4 -Identify top expressed genes" (Identify top expressed genes in spatial data and visualize them in a box plot.)
# INPUT OPTIONAL seurat_object.Robj: "Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL top_genes.pdf
# RUNTIME R-4.2.3-single-cell
# TOOLS_BIN ""

# 2022-07-22 IH

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(Matrix)

# Load the R-Seurat-object (called seurat_obj)
load("seurat_object.Robj")

# Open the pdf file for plotting
pdf(file = "top_genes.pdf", width = 13, height = 7)

# Identify top expressed genes
# code from  https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/seurat/seurat_07_spatial.html
C <- seurat_obj@assays$Spatial@counts
C@x <- C@x / rep.int(colSums(C), diff(C@p))
most_expressed <- order(Matrix::rowSums(C), decreasing = T)[20:1]
boxplot(as.matrix(t(C[most_expressed, ])),
    cex = 0.1, las = 1, xlab = "% total count per cell",
    col = (scales::hue_pal())(20)[20:1], horizontal = TRUE
)

# close the pdf
dev.off()

# EOF
