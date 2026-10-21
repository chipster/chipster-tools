# TOOL spatial-transcriptomics-seurat-find-markers-HD-v5.R: "Seurat v5 -Find marker genes for all clusters" (Given a Seurat object of Visium HD data with cluster information, this tool identifies marker genes for all clusters. It compares each cluster to all the other spots. The top 5 marker genes for each cluster are visualized as a heatmap.)
# INPUT seurat_obj_clustering.Robj: "Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL markers_Spatial.008um.pdf
# OUTPUT OPTIONAL markers_Spatial.016um.pdf
# OUTPUT OPTIONAL markers_Spatial.008um.tsv
# OUTPUT OPTIONAL markers_Spatial.016um.tsv
# PARAMETER assay: "Assay to use" TYPE [Spatial.008um: Spatial.008um, Spatial.016um: Spatial.016um] DEFAULT Spatial.008um
# PARAMETER banksy: "Do you want to use Banksy clusters" TYPE [TRUE: yes, FALSE: no] DEFAULT FALSE (Select yes if you identified spatially defined tissue domains with BANKSY.)
# PARAMETER OPTIONAL label.size: "determine the label size of the plots" TYPE INTEGER DEFAULT 3
# PARAMETER OPTIONAL width: "Width of the pdf" TYPE INTEGER DEFAULT 10
# PARAMETER OPTIONAL height: "Height of the pdf" TYPE INTEGER DEFAULT 10
# PARAMETER OPTIONAL min.pct: "Limit testing to genes which are expressed in at least this fraction of spots" TYPE DECIMAL DEFAULT 0.01 (Test only genes which are detected in at least this fraction of spots in either cluster. Withholding infrequently expressed genes will speed up testing.)
# PARAMETER OPTIONAL logfc.threshold: "Limit testing to genes which show at least this fold" TYPE DECIMAL DEFAULT 0.1 (Test only genes which show on average at least this log2 fold difference between the two groups of spots. Increasing the threshold speeds up testing, but can also miss weaker signals.)
# PARAMETER OPTIONAL test.use: "Test for differential expression" TYPE [wilcox: wilcox, MAST: MAST] DEFAULT wilcox
# PARAMETER OPTIONAL only.pos: "Report only positive marker genes" TYPE [TRUE: yes, FALSE: no] DEFAULT TRUE (If yes, only genes with positive log2 fold change are listed in the result file.)
# PARAMETER OPTIONAL p_val_adj: "Adjusted p-value cutoff" TYPE DECIMAL DEFAULT 0.1 (Only genes with an adjusted p-value smaller than this are kept in the result table)
# PARAMETER OPTIONAL min.pct: "Limit testing to genes which are expressed in at least this fraction of spots" TYPE DECIMAL DEFAULT 0.25 (Test only genes which are detected in at least this fraction of spots in either group of spots.)
# PARAMETER OPTIONAL logfc.threshold: "Limit testing to genes which show at least this fold change" TYPE DECIMAL DEFAULT 0.1 (Test only genes which show on average at least this log2 fold difference between the two groups of spots. Increasing the threshold speeds up testing, but can also miss weaker signals.)
# PARAMETER OPTIONAL label.size: "Label size of the plots" TYPE INTEGER DEFAULT 3
# PARAMETER OPTIONAL width: "Width of the plot" TYPE INTEGER DEFAULT 10
# PARAMETER OPTIONAL height: "Height of the plot" TYPE INTEGER DEFAULT 10
# RUNTIME R-4.5.1-visium-hd
# SLOTS 6
# TOOLS_BIN ""



# 2026-10-16 JV Add banksy clustering option

# This package was installed along the new docker image
#suppressMessages(install.packages("ape", repos = "https://cloud.r-project.org/"))

library("ape")
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyverse)
library("presto")

### PARAMETER reduction: "Reduction to use" TYPE [PCA] DEFAULT PCA
reduction <- as.character("PCA")
assay <- as.character(assay)
reduction <- tolower(as.character(reduction))
test.use <- as.character(test.use)
only.pos <- as.logical(only.pos)

load("seurat_obj_clustering.Robj")

#print(Reductions(seurat_obj))

# This is the one error still present, this "fixes" it
# for (g in Graphs(seurat_obj)) {
# seurat_obj[[g]] <- NULL
# }

# Tämä pitäisi saada toimimaan sitten. joskus oli "banksy.pca" reductioni.

if (banksy) {
  assay       <- "Spatial.008um"
  cluster_col <- "banksy_cluster"
  reduction   <- "BANKSY"          # check with Reductions(seurat_obj)
} else if (assay == "Spatial.008um") {
    cluster_col <- "seurat_cluster.projected"
    reduction   <- paste0("full.", reduction, ".sketch")
} else {
  cluster_col <- "seurat_cluster.016um"
  reduction   <- paste0(reduction, ".016um")
}



DefaultAssay(seurat_obj) <- assay
Idents(seurat_obj) <- cluster_col


# Downsampling cells may not needed here
#object_subset <- subset(seurat_obj, cells = Cells(seurat_obj[[assay]], downsample = 1000))


# Diet object for the analysis.
object_subset <- DietSeurat(object = seurat_obj, assays = assay, layers = c("counts", "data"), dimreducs = reduction, graphs = NULL)

object_subset <- subset(object_subset, cells = Cells(object_subset[[assay]]))

DefaultAssay(object_subset) <- assay
Idents(object_subset) <- cluster_col

object_subset <- BuildClusterTree(object_subset, assay = assay, reduction = reduction, reorder = TRUE, dims = 1:30)

print("Finding markers.")
markers <- FindAllMarkers(object_subset, assay = assay, only.pos = only.pos, test.use = test.use, min.pct = min.pct, logfc.threshold = logfc.threshold, print.bar = TRUE)
markers <- markers[markers$p_val_adj < p_val_adj,]

markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > logfc.threshold) %>%
  slice_head(n = 5) %>%
  ungroup() -> top5

object_subset <- ScaleData(object_subset, assay = assay, features = top5$gene)

pdf(file = paste0("markers_", assay, ".pdf"), width = width, height = height)

p <- DoHeatmap(object_subset, assay = assay, features = top5$gene, size = 2.5) + theme(axis.text = element_text(size = 5.5)) + NoLegend()

print(p)

dev.off()

save(markers, file = paste0("markers_", assay, ".tsv"))
write.table(as.matrix(markers), file = paste0("markers_", assay, ".tsv"), sep = "\t", row.names = T, col.names = T, quote = F)

# EOF