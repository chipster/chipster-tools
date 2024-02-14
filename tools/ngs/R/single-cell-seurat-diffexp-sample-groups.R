# TOOL single-cell-seurat-diffexp-sample-groups.R: "Seurat v4 -Find DE genes between sample groups" (This tool lists the differentially expressed genes between user defined sample groups for a user defined cluster. This tool can be used for Seurat objects with more than 2 samples in them.)
# INPUT OPTIONAL combined_seurat_obj.Robj: "Combined Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL de-list.tsv
# OUTPUT OPTIONAL de-list_{...}.tsv
# OUTPUT OPTIONAL expressionPlots.pdf
# PARAMETER OPTIONAL normalisation.method: "Normalisation method used previously" TYPE [LogNormalize:"Global scaling normalization", SCT:"SCTransform"] DEFAULT LogNormalize (Which normalisation method was used in preprocessing, Global scaling normalization \(default, NormalizeData function used\) or SCTransform.)
# PARAMETER samples1: "Name of the sample group to compare with" TYPE STRING DEFAULT "STIM" (Name of the sample group of which you want to identify the differentially expressed of.)
# PARAMETER samples2: "Name of the sample group to compare to" TYPE STRING DEFAULT "CTRL" (Name of the sample group which you want to identify the differentially expressed of.)
# PARAMETER cluster: "Name of the cluster" TYPE STRING DEFAULT 3 (Name of the cluster of which you want to identify the differentially expressed of. By default, the clusters are named with numbers starting from 0.)
# PARAMETER OPTIONAL only.positive: "Return only positive marker genes" TYPE [FALSE, TRUE] DEFAULT TRUE (Tool only returns positive markers as default. Change the parameter here if you want to also include the negative markers.)
# PARAMETER OPTIONAL logFC.de: "Fold change threshold for differentially expressed genes in log scale" TYPE DECIMAL FROM 0 TO 5 DEFAULT 0.25 (Genes with an average fold change smaller than this are not included in the analysis.)
# PARAMETER OPTIONAL pval.cutoff.de: "Adjusted p-value cutoff for differentially expressed genes" TYPE DECIMAL FROM 0 TO 1 DEFAULT 0.05 (Cutoff for the adjusted p-value of the DE genes: by default, adjusted p-values bigger than 0.05 are filtered out.)
# PARAMETER OPTIONAL minpct: "Limit testing for differentially expressed genes to genes which are expressed in at least this fraction of cells" TYPE DECIMAL DEFAULT 0.1 (Test only genes which are detected in at least this fraction of cells in either of two samples being compared in the cluster of question. Meant to speed up testing by leaving out genes that are very infrequently expressed.)
# RUNTIME R-4.2.3-single-cell
# SLOTS 2
# TOOLS_BIN ""

# 2018-16-05 ML
# 11.07.2019 ML Seurat v3
# 23.09.2019 EK Add only.pos = TRUE
# 30.10.2019 ML Add filtering parameters
# 02.12.2019 EK Change FC filter to use logfc.threshold prefiltering, add adjusted p-value filtering for DE genes
# 2021-10-04 ML Update to Seurat v4
# 2022-07-21 ML Tune for SCTransform data
# 2023-02-10 LG Add 2 slots
# 2023-04-03 ML Add parameters so that it is possible to print out all the genes + simplify code
# 2023-04-06 LG Remove 5 slots - Discrepancy in number of slots added v. removed
# 2023-07-11 ML Add 3 slots

library(Seurat)
library(dplyr)

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

# Add: Check that the sample group names are in the data!

# Set the identity as clusters:
sel.clust <- "seurat_clusters"
data.combined <- SetIdent(data.combined, value = sel.clust)
# table(data.combined@active.ident)

# Select all cells in cluster X (determined by user)
# NOTE: this could also be done using the subset.ident = "2" parameter in FindMarkers:
# # DGE_cell_selection <- FindMarkers(data.combined, ident.1 = samples1, ident.2 = samples2, group.by = 'type', subset.ident = cluster)
cell_selection <- subset(data.combined, cells = colnames(data.combined)[data.combined@meta.data[, sel.clust] == cluster])

# Set identity as "type" (this is where the sample group info is stored in setup tool)
cell_selection <- SetIdent(cell_selection, value = "type")


# Compute differentiall expression

if (normalisation.method == "SCT") {
  # Note: assay = "SCT" and recorrect_umi = FALSE
  DGE_cell_selection <- FindMarkers(cell_selection, assay = "SCT", ident.1 = samples1, ident.2 = samples2, group.by = "type", log2FC.threshold = logFC.de, min.pct = minpct, verbose = FALSE, return.thresh = pval.cutoff.de, recorrect_umi = FALSE, only.pos = only.positive) #, only.pos = only.positive) # min.diff.pct = 0.2, max.cells.per.ident = 50, test.use = "wilcox",

} else {
  DGE_cell_selection <- FindMarkers(cell_selection, assay = "RNA", ident.1 = samples1, ident.2 = samples2, group.by = "type", log2FC.threshold = logFC.de, min.pct = minpct,  verbose = FALSE, return.thresh = pval.cutoff.de, only.pos = only.positive) # min.diff.pct = 0.2, max.cells.per.ident = 50, test.use = "wilcox",
}
  
# Filter based on adj-p-val (no return.thresh parameter):
DGE_cell_selection_filtered <- DGE_cell_selection[DGE_cell_selection$p_val_adj<pval.cutoff.de, ]


# Add average expression to the table:
  if (normalisation.method == "SCT") {
    aver_expr <- AverageExpression(object = cell_selection, slot = "data", assay = "SCT")
  } else {
    aver_expr <- AverageExpression(object = cell_selection)
  }

  aver_expr_in_clusters <- aver_expr[[1]]
    
  # select the wanted columns (based on samples1.cluster and samples2.cluster ) and rows (DEGs):
  aver_expr_ident1 <- round(aver_expr_in_clusters[row.names(DGE_cell_selection_filtered),samples1 ], digits = 4)
  aver_expr_ident2 <- round(aver_expr_in_clusters[row.names(DGE_cell_selection_filtered),samples2 ], digits = 4)
  
  full_table <- cbind(DGE_cell_selection_filtered, aver_expr_ident1 , aver_expr_ident2)
  

# Comparison name for the output file:
comparison.name <- paste(samples1, "vs", samples2, "in_cluster", cluster, sep = "_")
name.for.output.file <- paste("de-list_", comparison.name, ".tsv", sep = "")

# Write to table
write.table(full_table, file = name.for.output.file, sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)


# Plots
pdf(file = "expressionPlots.pdf")

# Plot the expression across the “type”.

DGE_cell_selection %>%
  # group_by(cluster) %>%
  top_n(-5, p_val) -> top5_cell_selection

# VlnPlot(cell_selection, features = as.character(unique(top5_cell_selection$gene)), ncol = 3, group.by = "type", assay = "RNA", pt.size = 0.1)
VlnPlot(cell_selection, features = as.character(rownames(top5_cell_selection)), ncol = 3, group.by = "type", assay = "RNA", pt.size = 0.1)

# Plot these genes across all clusters, but split by “type”
# VlnPlot(data.combined, features = as.character(unique(top5_cell_selection$gene)), ncol = 3, split.by = "type", assay = "RNA", pt.size = 0)
VlnPlot(data.combined, features = as.character(rownames(top5_cell_selection)), ncol = 3, split.by = "type", assay = "RNA", pt.size = 0)

dev.off() # close the pdf

# Save the Robj for the next tool
# save(combined_seurat_obj, file="seurat_obj_combined.Robj")

## EOF
