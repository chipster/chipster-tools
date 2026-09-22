# TOOL single-cell-seurat-cluster-vis-v5.R: "Seurat v5 -Cluster visualisation" (Visualises clusters in tSNE and UMAP plots from a Seurat object that already contains cluster assignments.)
# INPUT seurat_obj.Robj: "Seurat object with cluster information" TYPE GENERIC
# OUTPUT OPTIONAL clusterPlot.pdf
# OUTPUT OPTIONAL log.txt
# PARAMETER OPTIONAL point.size: "Point size in tSNE and UMAP plots" TYPE DECIMAL DEFAULT 1 (Point size for the cluster plots.)
# PARAMETER OPTIONAL add.labels: "Add labels on top of clusters in plots" TYPE [TRUE: yes, FALSE: no] DEFAULT TRUE (Add cluster number on top of the cluster in UMAP and tSNE plots.)
# PARAMETER OPTIONAL reduction.to.plot: "Which dimensionality reductions to plot" TYPE [umap_and_tsne: "UMAP and tSNE", umap_only: "UMAP only", tsne_only: "tSNE only"] DEFAULT umap_and_tsne (Select which dimensionality reduction plots to produce. If the reduction was not computed in the clustering step, that plot will be skipped.)
# RUNTIME R-4.3.2-single-cell
# TOOLS_BIN ""


# 2026-09-18 ML Split off from single-cell-seurat-clustering-v5.R for standalone use

library(Seurat)
library(gplots)
options(Seurat.object.assay.version = "v5")

# Load the Seurat object (may be named seurat_obj or data.combined)
load("seurat_obj.Robj")

if (exists("data.combined")) {
  seurat_obj <- data.combined
}

# Sanity check: clusters must exist
if (is.null(Idents(seurat_obj)) || length(levels(Idents(seurat_obj))) == 0) {
  stop("CHIPSTER-NOTE: The Seurat object does not appear to contain cluster assignments. Please run the clustering tool first.")
}

# Check which reductions are available
available_reductions <- names(seurat_obj@reductions)
has_umap <- "umap" %in% available_reductions
has_tsne <- "tsne" %in% available_reductions

if (!has_umap && !has_tsne) {
  stop("CHIPSTER-NOTE: The Seurat object contains neither a UMAP nor a tSNE reduction. Please run the clustering tool (which computes both) before using this visualisation tool.")
}

# Determine what to plot based on parameter and availability
plot_umap <- (reduction.to.plot %in% c("umap_and_tsne", "umap_only")) && has_umap
plot_tsne <- (reduction.to.plot %in% c("umap_and_tsne", "tsne_only")) && has_tsne

if (!plot_umap && !plot_tsne) {
  stop(paste0("CHIPSTER-NOTE: None of the requested reductions are available in this object. ",
              "Available reductions: ", paste(available_reductions, collapse = ", "), "."))
}

# Warn if a requested reduction is missing (but still produce the available one)
if ((reduction.to.plot %in% c("umap_and_tsne", "umap_only")) && !has_umap) {
  message("Warning: UMAP reduction not found in the Seurat object — skipping UMAP plot.")
}
if ((reduction.to.plot %in% c("umap_and_tsne", "tsne_only")) && !has_tsne) {
  message("Warning: tSNE reduction not found in the Seurat object — skipping tSNE plot.")
}

# Open PDF
pdf(file = "clusterPlot.pdf")

if (plot_umap) {
  DimPlot(seurat_obj, reduction = "umap", pt.size = point.size, label = add.labels)
}

if (plot_tsne) {
  DimPlot(seurat_obj, reduction = "tsne", pt.size = point.size, label = add.labels)
}

# Number of cells in each cluster
cell_counts <- table(Idents(seurat_obj), seurat_obj$orig.ident)
textplot(cell_counts, halign = "center", valign = "center", cex = 1)
title(paste("Total number of cells:", length(colnames(x = seurat_obj)),
            "\nNumber of cells in each cluster:"))

dev.off()

# EOF