# TOOL spatial-transcriptomics-seurat-find-markers-HD-v5.R: "Seurat v5 HD -Find markers" (This tool identifies marker genes.)
# INPUT seurat_obj_clustering.Robj: "Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL markers_Spatial.008um.pdf
# OUTPUT OPTIONAL markers_Spatial.016um.pdf
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

suppressMessages(install.packages("ape", repos = "https://cloud.r-project.org/"))

library("ape")
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyverse)

load("seurat_obj_clustering.Robj")

if ("sketch" %in% Assays(seurat_obj)) {
  seurat_obj[["sketch"]] <- NULL
}

for (g in Graphs(seurat_obj)) {
seurat_obj[[g]] <- NULL
}

DefaultAssay(seurat_obj) <- "Spatial.008um"
reduction <- Reductions(seurat_obj)[6] #Full umap sketch

print("colnames of meta.data")
print(colnames(seurat_obj@meta.data))


head(Cells(seurat_obj[["Spatial.008um"]]))

head(Reductions(seurat_obj))


#Idents(object) <- "seurat_cluster.projected"
#object_subset <- subset(seurat_obj, cells = Cells(seurat_obj[[assay]]), downsample=1000)

object_subset <- seurat_obj


#print(reduction)



## Order clusters by similarity
#DefaultAssay(object_subset) <- assay
Idents(object_subset) <- "seurat_cluster.projected"

head(Idents(object_subset))

object_subset <- BuildClusterTree(object_subset, assay = "Spatial.008um", reduction = reduction, reorder = T)

FindMarkers(object_subset, assay = "Spatial.008um", ident.1 = 0)

#markers <- FindAllMarkers(object_subset, assay = "Spatial.008um")
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup() -> top5

object_subset <- ScaleData(object_subset, assay = "Spatial.008um", features = top5$gene)

pdf(file = paste0("markers_", "Spatial.008um", ".pdf"), width = width, height = height)

p <- DoHeatmap(object_subset, assay = "Spatial.008um", features = top5$gene, size = 2.5) + theme(axis.text = element_text(size = 5.5)) + NoLegend()

print(p)


# EOF