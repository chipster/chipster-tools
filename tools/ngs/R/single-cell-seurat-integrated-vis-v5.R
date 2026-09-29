# TOOL single-cell-seurat-integrated-vis-v5.R: "Seurat v5 -Integrated analysis visualisation" (Visualises clusters from an integrated multi-sample Seurat object. Produces UMAP plots coloured by sample and cluster identity, and a cell-count summary table. To be used after the Seurat v5 -Integrated analysis tool.)
# INPUT seurat_obj_integrated.Robj: "Integrated Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL integrated_plot.pdf
# OUTPUT OPTIONAL aver_expr_in_clusters.tsv
# OUTPUT OPTIONAL log.txt
# PARAMETER OPTIONAL point.size: "Point size in feature plots" TYPE DECIMAL DEFAULT 1 (Point size for the feature plots.)
# PARAMETER OPTIONAL width: "Page width" TYPE INTEGER DEFAULT 24 (Width of the output plots in inches.)
# PARAMETER OPTIONAL height: "Height of the output plots" TYPE INTEGER DEFAULT 12 (Height of the output plots in inches.)
# PARAMETER OPTIONAL add.labels: "Add labels on top of clusters in plots" TYPE [TRUE: yes, FALSE: no] DEFAULT TRUE (Add cluster number on top of each cluster in the UMAP plots.)
# PARAMETER OPTIONAL output_aver_expr: "Give a list of average expression in each cluster" TYPE [T: yes, F: no] DEFAULT F (Returns an expression table for an average single cell in each cluster.)
# PARAMETER OPTIONAL normalisation.method: "Normalisation method used previously" TYPE [LogNormalize:"Global scaling normalization", SCT:"SCTransform"] DEFAULT LogNormalize (Which normalisation method was used in preprocessing. This is needed only when requesting the average expression table.)
# RUNTIME R-4.3.2-single-cell
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
has_annotations       <- "seurat_annotations" %in% colnames(data.combined@meta.data)

# --- Layout and scaling ---
# n_cols_split: how many split-by-sample panels fit side-by-side with ~square aspect ratio.
# Minimum 2 to avoid single-panel pages on a square page.
samples      <- levels(factor(data.combined$stim))
n_cols_split <- max(2, floor(width / height))
max_per_page <- n_cols_split   # one row of panels per split page

# Base text size on the smaller of panel width vs height.
# For 2-panel pages each panel is width/2, so that is the constraining dimension.
# This prevents long titles from overflowing narrow panels.
panel_size <- min(width / 2, height)
text_size  <- round(11 * panel_size / 6)

scale_text <- function(sz) {
  theme(
    plot.title   = element_text(size = sz),
    axis.title   = element_text(size = sz),
    axis.text    = element_text(size = sz * 0.8),
    legend.text  = element_text(size = sz * 0.8),
    legend.title = element_text(size = sz)
  )
}

# Scale textplot cex with page height so the table fills the page.
cex_table <- max(1.0, 1.2 * height / 12)

pdf(file = "integrated_plot.pdf", width = width, height = height)


# --- Page 1: Unintegrated UMAP, colour by sample and by cluster ---
if (has_umap_unintegrated) {
  p1 <- DimPlot(data.combined, reduction = "umap.unintegrated",
                group.by = "stim", pt.size = point.size) +
        ggtitle("UMAP unintegrated, color by sample") +
        scale_text(text_size)
  p2 <- DimPlot(data.combined, reduction = "umap.unintegrated",
                group.by = "seurat_clusters", pt.size = point.size, label = add.labels) +
        ggtitle("UMAP unintegrated, color by cluster") +
        scale_text(text_size)
  print(plot_grid(p1, p2))
}

# --- Page 2: Integrated UMAP, colour by sample and by cluster ---
p3 <- DimPlot(data.combined, reduction = "umap",
              group.by = "stim", pt.size = point.size) +
      ggtitle("Chosen reduction method, integrated, color by sample") +
      scale_text(text_size)

if (has_annotations) {
  p4 <- DimPlot(data.combined, reduction = "umap",
                group.by = "seurat_annotations", pt.size = point.size, label = add.labels) +
        ggtitle("Chosen reduction method, integrated, color by cluster") +
        scale_text(text_size)
} else {
  p4 <- DimPlot(data.combined, reduction = "umap",
                pt.size = point.size, label = add.labels) +
        ggtitle("Chosen reduction method, integrated, color by cluster") +
        scale_text(text_size)
}
print(plot_grid(p3, p4))

# --- Page 3+: Integrated UMAP split by sample (paginated if many samples) ---
if (length(samples) <= max_per_page) {
  p5 <- DimPlot(data.combined, reduction = "umap",
                split.by = "stim", ncol = n_cols_split,
                pt.size = point.size, label = add.labels) +
        ggtitle("Chosen reduction method, integrated, samples") +
        scale_text(text_size)
  print(p5)
} else {
  batches <- split(samples, ceiling(seq_along(samples) / max_per_page))
  for (i in seq_along(batches)) {
    obj_sub <- subset(data.combined, stim %in% batches[[i]])
    # Keep original cluster levels so colours stay consistent across pages
    Idents(obj_sub) <- factor(Idents(obj_sub), levels = levels(Idents(data.combined)))
    p5 <- DimPlot(obj_sub, reduction = "umap",
                  split.by = "stim", ncol = n_cols_split,
                  pt.size = point.size, label = add.labels) +
          ggtitle(sprintf("Chosen reduction method, integrated, samples (%d/%d)",
                          i, length(batches))) +
          scale_text(text_size)
    print(p5)
  }
}

 # --- Cell count table ---
cell_counts <- table(Idents(data.combined), data.combined$stim)
sums <- colSums(cell_counts)
cell_counts <- rbind(cell_counts, sums)

textplot(cell_counts, halign = "center", valign = "center", cex = 2)
title(paste("Total number of cells:", length(colnames(x = data.combined)), "\nNumber of cells in each cluster:"))

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