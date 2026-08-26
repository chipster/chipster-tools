# TOOL spatial-transcriptomics-seurat-find-markers-HD-v5.R: "Seurat v5 HD -Find markers" (This tool identifies marker genes.)
# INPUT seurat_obj_clustering.Robj: "Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL markers_Spatial.008um.pdf
# OUTPUT OPTIONAL markers_Spatial.016um.pdf
# PARAMETER assay: "Assay to use" TYPE [Spatial.008um: Spatial.008um, Spatial.016um: Spatial.016um] DEFAULT Spatial.008um
# PARAMETER OPTIONAL label.size: "determine the label size of the plots" TYPE INTEGER DEFAULT 3
# PARAMETER OPTIONAL width: "Width of the pdf" TYPE INTEGER DEFAULT 10
# PARAMETER OPTIONAL height: "Height of the pdf" TYPE INTEGER DEFAULT 10
# PARAMETER OPTIONAL min.pct: "Limit testing to genes which are expressed in at least this fraction of spots" TYPE DECIMAL DEFAULT 0.01 (Test only genes which are detected in at least this fraction of spots in either cluster. Withholding infrequently expressed genes will speed up testing.)
# PARAMETER OPTIONAL logfc.threshold: "Limit testing to genes which show at least this fold" TYPE DECIMAL DEFAULT 0.1 (Test only genes which show on average at least this log2 fold difference between the two groups of spots. Increasing the threshold speeds up testing, but can also miss weaker signals.)
# PARAMETER OPTIONAL test.use: "Test for differential expression" TYPE [wilcox: wilcox, MAST: MAST] DEFAULT wilcox
# PARAMETER OPTIONAL only.pos: "Report only positive marker genes" TYPE [FALSE, TRUE] DEFAULT FALSE (By default, this tool)
# RUNTIME R-4.5.1-seurat5
# SLOTS 10
# TOOLS_BIN ""


# This package was installed along the new docker image
#suppressMessages(install.packages("ape", repos = "https://cloud.r-project.org/"))

assay <- as.character(assay)

library("ape")
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyverse)
library("presto")

load("seurat_obj_clustering.Robj")

print("Nimet!!")
#print((Reductions(seurat_obj)))

print("Colnames")
#print(colnames(seurat_obj@meta.data))


# Maybe a useless step???
# if ("sketch" %in% Assays(seurat_obj)) {
#   seurat_obj[["sketch"]] <- NULL
# }


# This is the one error still present, this "fixes" it
for (g in Graphs(seurat_obj)) {
seurat_obj[[g]] <- NULL
}

# Crete downsampled object to make visualization either
if (assay == "Spatial.008um") {
  cluster_col <- "seurat_cluster.008um"
  reduction   <- "pca.008um"
} else {
  cluster_col <- "seurat_cluster.016um"
  reduction   <- "pca.016um"
}

DefaultAssay(seurat_obj) <- assay
Idents(seurat_obj) <- cluster_col

object_subset <- seurat_obj
DefaultAssay(object_subset) <- assay
Idents(object_subset) <- cluster_col

object_subset <- BuildClusterTree(object_subset, assay = assay, reduction = reduction, reorder = TRUE)

print("Finding markers...")
markers <- FindAllMarkers(object_subset, assay = assay, only.pos = TRUE, print.bar = TRUE)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0.5) %>%
  slice_head(n = 5) %>%
  ungroup() -> top5

object_subset <- ScaleData(object_subset, assay = assay, features = top5$gene)

pdf(file = paste0("markers_", assay, ".pdf"), width = width, height = height)

p <- DoHeatmap(object_subset, assay = assay, features = top5$gene, size = 2.5) + theme(axis.text = element_text(size = 5.5)) + NoLegend()

print(p)

dev.off()


# EOF
