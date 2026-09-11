# TOOL spatial-transcriptomics-seurat-diffexp-chosen-clusters-HD-v5.R: "Seurat v5 HD -Identify spatially variable genes based on clusters" (This tool identifies differentially expressed genes between two user defined clusters and visualizes these genes on top of the tissue image.)
# INPUT seurat_obj_clustering.Robj: "Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL Markerplot.pdf
# OUTPUT OPTIONAL spatially_variable_genes.tsv
# PARAMETER assay: "Assay to use" TYPE [Spatial.008um: Spatial.008um, Spatial.016um: Spatial.016um] DEFAULT Spatial.008um
# PARAMETER OPTIONAL sketch: "Was sketch clustering used for 8 bin assay" TYPE [TRUE: yes, FALSE: no] DEFAULT FALSE ()
# PARAMETER OPTIONAL cluster1: "First cluster" TYPE INTEGER DEFAULT 1 (Cluster you want to identify the differentially expressed for.)
# PARAMETER OPTIONAL cluster2: "Second cluster" TYPE INTEGER DEFAULT 2 (A second cluster for comparison.)
# PARAMETER OPTIONAL min_pct: "Limit testing to genes which are expressed in at least this fraction of spots" TYPE DECIMAL DEFAULT 0.01 (Test only genes which are detected in at least this fraction of spots in either cluster. Withholding infrequently expressed genes will speed up testing.)
# PARAMETER OPTIONAL logfc_threshold: "Limit testing to genes which show at least this fold difference" TYPE DECIMAL DEFAULT 0.1 (Test only genes which show on average at least this log2 fold difference between the two groups of spots. Increasing the threshold speeds up testing, but can also miss weaker signals.)
# PARAMETER OPTIONAL test: "Test for differential expression" TYPE [wilcox: wilcox, MAST: MAST] DEFAULT wilcox
# PARAMETER OPTIONAL only_pos: "Report only positive marker genes" TYPE [FALSE, TRUE] DEFAULT FALSE (By default, this tool lists all markers. When this parameter is set to TRUE, only genes with positive log2 fold change are listed in the result file.)
# PARAMETER OPTIONAL no_of_feats: "Number of spatially variable genes to visualize" TYPE INTEGER DEFAULT 3 (Choose the number of highest variable genes to visualize.)
# PARAMETER OPTIONAL color.scale: "Determine color scale based on all genes" TYPE [all:yes, feature:no] DEFAULT feature (Determine whether the color scale is based on all genes or individual genes. By default, the color scale is determined for each gene individually and may differ between genes.)
# PARAMETER OPTIONAL width: "Width of the output pdf" TYPE DECIMAL DEFAULT 10 (Width of the output pdf in inches)
# PARAMETER OPTIONAL height: "Height of the output pdf" TYPE DECIMAL DEFAULT 10 (Height of the output pdf in inches)
# RUNTIME R-4.5.1-seurat5
# SLOTS 4
# TOOLS_BIN ""



# 2022-07-29 IH
# 2022-10-20 ML Add output for spatially_variable_genes.tsv
# 2023-02-23 LG Add parameters min.pct, logfc.threshold, and only.pos
# 2024-03-21 EP Update to Seurat v5
# 2026-09-10 JV Update to Visium HD

assay <- as.character(assay)
sketch <- as.logical(sketch)
cluster1 <- as.numeric(cluster1)
cluster2 <- as.numeric(cluster2)
min_pct <- as.numeric(min_pct)
logfc_threshold <- as.numeric(logfc_threshold)
test <- as.character(test)
only_pos <- as.logical(only_pos)


library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(Biobase)

print(package.version("Seurat"))
documentVersion("Seurat", package.version("Seurat"))

# Load the R-Seurat-object (called seurat_obj)
load("seurat_obj_clustering.Robj")



if(assay == "Spatial.008um") {
  DefaultAssay(seurat_obj) <- "Spatial.008um"

  if(sketch) {
    Idents(seurat_obj) <- seurat_obj$seurat_cluster.projected
  } else {
    Idents(seurat_obj) <- seurat_obj$seurat_cluster.008um
  }
} else {
  DefaultAssay(seurat_obj) <- "Spatial.016um"

  Idents(seurat_obj) <- seurat_obj$seurat_cluster.016um
}

#print(Idents(seurat_obj))


# SCTransform not good for Visium HD
# Assay Spatial 8bin or 16bin shall be used instead
#seurat_obj <- PrepSCTFindMarkers(object = seurat_obj, assay = "SCT")

print(DefaultAssay(seurat_obj))
# Differential expression
de_markers <- FindMarkers(seurat_obj, assay = assay, ident.1 = cluster1, ident.2 = cluster2, test.use = test, logfc.threshold = logfc_threshold, min.pct = min_pct, only.pos = only_pos)

# Print out markers into a table:
# name.for.file <- paste("spatially_variable_genes_cluster", cluster1, "_vs_cluster", cluster2, ".tsv", sep="")
write.table(as.matrix(de_markers), file = "spatially_variable_genes.tsv", sep = "\t", row.names = T, col.names = T, quote = F)

# Open the pdf file for plotting
pdf(file = "Markerplot.pdf", , width = width, height = height)

SpatialFeaturePlot(object = seurat_obj, features = rownames(de_markers)[1:no_of_feats], keep.scale = color.scale, alpha = c(0.1, 1))

dev.off() # close the pdf

# EOF