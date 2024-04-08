# TOOL single-cell-seurat-diffexp-pseudobulk-v5.R: "Seurat v5 -Find DE genes between sample groups, pseudobulk" (This tool lists the differentially expressed genes between user defined sample groups for a user defined cluster. This tool can be used for Seurat objects with more than 2 samples in them.)
# INPUT OPTIONAL combined_seurat_obj.Robj: "Combined Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL de-list.tsv
# OUTPUT OPTIONAL pseudobulk_de-list_{...}.tsv
# OUTPUT OPTIONAL expressionPlots.pdf
# PARAMETER OPTIONAL normalisation.method: "Normalisation method used previously" TYPE [LogNormalize:"Global scaling normalization", SCT:"SCTransform"] DEFAULT LogNormalize (Which normalisation method was used in preprocessing, Global scaling normalization \(default, NormalizeData function used\) or SCTransform.)
# PARAMETER samples1: "Name of the sample group to compare with" TYPE STRING DEFAULT "STIM" (Name of the sample group of which you want to identify the differentially expressed of.)
# PARAMETER samples2: "Name of the sample group to compare to" TYPE STRING DEFAULT "CTRL" (Name of the sample group of which you want to identify the differentially expressed of.)
# PARAMETER cluster: "Name of the cluster" TYPE STRING DEFAULT 3 (Name of the cluster of which you want to identify the differentially expressed of. By default, the clusters are named with numbers starting from 0.)
# PARAMETER OPTIONAL only.positive: "Return only positive marker genes" TYPE [FALSE, TRUE] DEFAULT TRUE (Tool only returns positive markers as default. Change the parameter here if you want to also include the negative markers.)
# PARAMETER OPTIONAL logFC.de: "Fold change threshold for differentially expressed genes in log2 scale" TYPE DECIMAL FROM 0 TO 5 DEFAULT 0.25 (Genes with an average fold change smaller than this are not included in the analysis.)
# PARAMETER OPTIONAL pval.cutoff.de: "Adjusted p-value cutoff for differentially expressed genes" TYPE DECIMAL FROM 0 TO 1 DEFAULT 0.05 (Cutoff for the adjusted p-value of the DE genes: by default, adjusted p-values bigger than 0.05 are filtered out.)
# PARAMETER OPTIONAL minpct: "Limit testing for differentially expressed genes to genes which are expressed in at least this fraction of cells" TYPE DECIMAL DEFAULT 0.1 (Test only genes which are detected in at least this fraction of cells in either of two samples being compared in the cluster of question. Meant to speed up testing by leaving out genes that are very infrequently expressed.)
# RUNTIME R-4.3.2-single-cell
# SLOTS 2
# TOOLS_BIN ""

# 2024-03-07 ML

# Based on: https://satijalab.org/seurat/articles/parsebio_sketch_integration

library(Seurat)
library(dplyr)
library(ggrepel)
library(ggplot2)
library(cowplot)
options(Seurat.object.assay.version = "v5")

# Load the R-Seurat-objects
load("combined_seurat_obj.Robj")

if (exists("seurat_obj")) {
  data.combined <- seurat_obj
}

DefaultAssay(data.combined) <- "RNA" # this is very crucial.

if (normalisation.method == "SCT") {
  # When SCTransform was used to normalise the data, do a prep step:
  data.combined <- PrepSCTFindMarkers(data.combined)
}

# Set the identity as clusters:
sel.clust <- "seurat_clusters"
data.combined <- SetIdent(data.combined, value = sel.clust)
# table(data.combined@active.ident)


# Aggregate 
bulk <- AggregateExpression(data.combined,  return.seurat = T, assays = "RNA",  group.by = c("seurat_clusters", "stim", "type")) 
# Note: "First group.by variable `seurat_clusters` starts with a number, appending `g` to ensure valid variable names"
tail(Cells(bulk))
# [1] "g13_normal5_NORMAL"  "g14_covid1_COVID"    "g14_covid17_COVID"   "g14_normal13_NORMAL"

# Subset based on cluster number
cluster.bulk <- subset(bulk, seurat_clusters == paste("g",cluster, sep="")) 
# Note, "g" added by AggregateExpression

# "type" contains the sample group information from the setup tool
Idents(cluster.bulk) <- "type"


# Compute differential expression
if (normalisation.method == "SCT") {
  DGE_cell_selection <- FindMarkers(cell_selection, assay = "SCT", ident.1 = samples1, ident.2 = samples2, group.by = "type", log2FC.threshold = logFC.de, min.pct = minpct, verbose = FALSE, return.thresh = pval.cutoff.de, recorrect_umi = FALSE, only.pos = only.positive) #, only.pos = only.positive) # min.diff.pct = 0.2, max.cells.per.ident = 50, test.use = "wilcox",
  pseudobulk_markers <- FindMarkers(cluster.bulk, assay = "SCT", ident.1 = samples1, ident.2 = samples2, slot = "counts", test.use = "DESeq2", verbose = F)
  # kesken!
} else {
  pseudobulk_markers <- FindMarkers(cluster.bulk, ident.1 = samples1, ident.2 = samples2, slot = "counts", test.use = "DESeq2", verbose = F)
}
 
# Filter based on adj-p-val (no return.thresh parameter):
pseudobulk_markers <- pseudobulk_markers[pseudobulk_markers$p_val_adj<pval.cutoff.de, ]

# Plot
pdf(file = "expressionPlots.pdf")

pseudobulk_markers$gene <- rownames(pseudobulk_markers)
ggplot(pseudobulk_markers, aes(avg_log2FC, -log10(p_val))) + geom_point(size = 0.5, alpha = 0.5) + theme_bw() +
  ylab("-log10(unadjusted p-value)") + geom_text_repel(aes(label = ifelse(p_val_adj < 0.01, gene,
                                                                          "")), colour = "red", size = 3)

pseudobulk_markers %>%
  # group_by(cluster) %>%
  top_n(-6, p_val) -> top5_cell_selection

# VlnPlot(cell_selection, features = as.character(unique(top5_cell_selection$gene)), ncol = 3, group.by = "type", assay = "RNA", pt.size = 0.1)
VlnPlot(cluster.bulk, features = as.character(rownames(top5_cell_selection)), ncol = 3, group.by = "type", assay = "RNA", pt.size = 0.1) + ggtitle("pseudobulk top genes")

# Plot these genes across all clusters, but split by “type”
# VlnPlot(data.combined, features = as.character(unique(top5_cell_selection$gene)), ncol = 3, split.by = "type", assay = "RNA", pt.size = 0)
VlnPlot(data.combined, features = as.character(rownames(top5_cell_selection)), ncol = 3, split.by = "type", assay = "RNA", pt.size = 0) + ggtitle("all cells top genes, clusters")

Idents(data.combined) <- "stim"
#VlnPlot(data.combined, features = as.character(rownames(top5_cell_selection)), ncol = 3, split.by = "stim", assay = "RNA", pt.size = 0)
VlnPlot(data.combined, features = as.character(rownames(top5_cell_selection)), ncol = 3, split.by = "type", assay = "RNA", pt.size = 0)  + ggtitle("all cells top genes, samples")


# According: https://nbisweden.github.io/workshop-scRNAseq/labs/seurat/seurat_05_dge.html#pseudobulk
pseudobulk_markers %>%
  # group_by(cluster) %>%
  top_n(-10, p_val) -> top10_cell_selection

DotPlot(cluster.bulk,
        features = top10_cell_selection$gene, group.by = "orig.ident",
        assay = "RNA"
) + coord_flip() + ggtitle("EdgeR pseudobulk") + RotatedAxis()

DotPlot(data.combined,
        features = top10_cell_selection$gene, group.by = "stim",
        assay = "RNA"
) + coord_flip() + ggtitle("EdgeR all cells") + RotatedAxis()

dev.off() # close the pdf



# # Add average expression to the table:
#   if (normalisation.method == "SCT") {
#     aver_expr <- AverageExpression(object = cell_selection, slot = "data", assay = "SCT")
#   } else {
#     aver_expr <- AverageExpression(object = cell_selection)
#   }

#   aver_expr_in_clusters <- aver_expr[[1]]
    
#   # select the wanted columns (based on samples1.cluster and samples2.cluster ) and rows (DEGs):
#   aver_expr_ident1 <- round(aver_expr_in_clusters[row.names(DGE_cell_selection_filtered),samples1 ], digits = 4)
#   aver_expr_ident2 <- round(aver_expr_in_clusters[row.names(DGE_cell_selection_filtered),samples2 ], digits = 4)
  
#   full_table <- cbind(DGE_cell_selection_filtered, aver_expr_ident1 , aver_expr_ident2)
  

# Comparison name for the output file:
comparison.name <- paste(samples1, "vs", samples2, "in_cluster", cluster, sep = "_")
name.for.output.file <- paste("pseudobulk_de-list_", comparison.name, ".tsv", sep = "")

# Write to table
write.table(pseudobulk_markers, file = name.for.output.file, sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)


# # Plots
# pdf(file = "expressionPlots.pdf")

# # Plot the expression across the “type”.

# DGE_cell_selection %>%
#   # group_by(cluster) %>%
#   top_n(-5, p_val) -> top5_cell_selection

# # VlnPlot(cell_selection, features = as.character(unique(top5_cell_selection$gene)), ncol = 3, group.by = "type", assay = "RNA", pt.size = 0.1)
# VlnPlot(cell_selection, features = as.character(rownames(top5_cell_selection)), ncol = 3, group.by = "type", assay = "RNA", pt.size = 0.1)

# # Plot these genes across all clusters, but split by “type”
# # VlnPlot(data.combined, features = as.character(unique(top5_cell_selection$gene)), ncol = 3, split.by = "type", assay = "RNA", pt.size = 0)
# VlnPlot(data.combined, features = as.character(rownames(top5_cell_selection)), ncol = 3, split.by = "type", assay = "RNA", pt.size = 0)

# dev.off() # close the pdf

# Save the Robj for the next tool
# save(combined_seurat_obj, file="seurat_obj_combined.Robj")

## EOF
