# TOOL single-cell-seurat-integrated-vis-v5.R: "Seurat v5 -Integrated analysis visualisation" (Visualises clusters from an integrated multi-sample Seurat object. Produces UMAP plots coloured by sample and cluster identity, and a cell-count summary table. To be used after the Seurat v5 -Integrated analysis tool.)
# INPUT seurat_obj_integrated.Robj: "Integrated Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL integrated_plot.pdf
# OUTPUT OPTIONAL aver_expr_in_clusters.tsv
# OUTPUT OPTIONAL log.txt
# PARAMETER OPTIONAL point.size: "Point size in plots" TYPE DECIMAL DEFAULT 0.5 (Point size for the dimensionality reduction plots.)
# PARAMETER OPTIONAL add.labels: "Add labels on top of clusters in plots" TYPE [TRUE: yes, FALSE: no] DEFAULT TRUE (Add cluster number on top of each cluster in the UMAP plots.)
# PARAMETER OPTIONAL output_aver_expr: "Give a list of average expression in each cluster" TYPE [T: yes, F: no] DEFAULT F (Returns an expression table for an average single cell in each cluster.)
# PARAMETER OPTIONAL normalisation.method: "Normalisation method used previously" TYPE [LogNormalize:"Global scaling normalization", SCT:"SCTransform"] DEFAULT LogNormalize (Which normalisation method was used in preprocessing. This is needed only when requesting the average expression table.)
# RUNTIME R-4.3.1-single-cell
# TOOLS_BIN ""

# 2026-09-18 ML Split off from single-cell-seurat-integrated-analysis-v5.R for standalone use

library(Seurat)
library(gplots)
library(ggplot2)
require(cowplot)
options(Seurat.object.assay.version = "v5")

# Load the integrated Seurat object (named data.combined)
load("seurat_obj_integrated.Robj")

# Sanity checks
if (!exists("data.combined")) {
  stop("CHIPSTER-NOTE: The expected object 'data.combined' was not found in the Robj file. Please make sure you are providing the output of the Seurat v5 -Integrated analysis tool.")
}

if (is.null(Idents(data.combined)) || length(levels(Idents(data.combined))) == 0) {
  stop("CHIPSTER-NOTE: The Seurat object does not contain cluster assignments. Please run the Seurat v5 -Integrated analysis tool first.")
}

available_reductions <- names(data.combined@reductions)

if (!"umap" %in% available_reductions) {
  stop(paste0("CHIPSTER-NOTE: No UMAP reduction found in the object. Available reductions: ",
              paste(available_reductions, collapse = ", "), ". Please run the Seurat v5 -Integrated analysis tool first."))
}

# Check whether the unintegrated UMAP and seurat_annotations are available
# (they may not be present in all versions of the pipeline)
has_umap_unintegrated <- "umap.unintegrated" %in% available_reductions
has_annotations <- "seurat_annotations" %in% colnames(data.combined@meta.data)

# Open PDF (wide format to fit side-by-side panels)
pdf(file = "integrated_plot.pdf", width = 13, height = 7)

# --- Plot 1 & 2: unintegrated vs integrated UMAP (side by side) ---
if (has_umap_unintegrated) {
  p1 <- DimPlot(data.combined, reduction = "umap.unintegrated",
                group.by = c("stim", "seurat_clusters"),
                pt.size = point.size)
} else {
  message("Warning: 'umap.unintegrated' reduction not found — substituting integrated UMAP coloured by sample.")
  p1 <- DimPlot(data.combined, reduction = "umap",
                group.by = "stim",
                pt.size = point.size)
}

if (has_annotations) {
  p2 <- DimPlot(data.combined, reduction = "umap",
                group.by = c("stim", "seurat_annotations"),
                pt.size = point.size, label = add.labels)
} else {
  message("Warning: 'seurat_annotations' metadata column not found — colouring by cluster identity instead.")
  p2 <- DimPlot(data.combined, reduction = "umap",
                group.by = "seurat_clusters",
                pt.size = point.size, label = add.labels)
}

# --- Plot 3: integrated UMAP split by sample ---
p3 <- DimPlot(data.combined, reduction = "umap",
              split.by = "stim",
              pt.size = point.size, label = add.labels)

print(plot_grid(p1, p2))
print(p3)

# --- Cell count table ---
cell_counts <- table(Idents(data.combined), data.combined$stim)
sums <- colSums(cell_counts)
cell_counts <- rbind(cell_counts, sums)

textplot(cell_counts, halign = "center", valign = "center", cex = 1)
title(paste("Total number of cells:", length(colnames(x = data.combined)),
            "\nNumber of cells in each cluster:"))

dev.off()

# Average expression table
# If requested, return expression for an 'average' single cell in each cluster.
if (output_aver_expr == "T") {
  if (normalisation.method == "SCT") {
    aver_expr <- AverageExpression(object = data.combined, slot = "data", assay = "SCT")
  } else {
    aver_expr <- AverageExpression(object = data.combined)
  }
  aver_expr_in_clusters <- aver_expr[[1]]
  write.table(aver_expr_in_clusters, file = "aver_expr_in_clusters.tsv", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
}

# EOF