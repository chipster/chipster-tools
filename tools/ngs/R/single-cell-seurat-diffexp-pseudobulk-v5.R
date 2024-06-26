# TOOL single-cell-seurat-diffexp-pseudobulk-v5.R: "Seurat v5 -Find DE genes between sample groups, pseudobulk" (This tool lists the differentially expressed genes between user defined sample groups for a user defined cluster. This tool can be used for Seurat objects with more than 2 samples in them.)
# INPUT OPTIONAL combined_seurat_obj.Robj: "Combined Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL de-list.tsv
# OUTPUT OPTIONAL pseudobulk_de-list_{...}.tsv
# OUTPUT OPTIONAL expressionPlots.pdf
# OUTPUT OPTIONAL pseudobulk_markers_NAs.tsv
# OUTPUT OPTIONAL dot_plot.Robj
# OUTPUT OPTIONAL data_combined_pseudobulk.Robj
# PARAMETER OPTIONAL normalisation.method: "Normalisation method used previously" TYPE [LogNormalize:"Global scaling normalization", SCT:"SCTransform"] DEFAULT LogNormalize (Which normalisation method was used in preprocessing, Global scaling normalization \(default, NormalizeData function used\) or SCTransform.)
# PARAMETER samples1: "Name of the sample group to compare with" TYPE STRING DEFAULT "STIM" (Name of the sample group of which you want to identify the differentially expressed genes of.)
# PARAMETER samples2: "Name of the sample group to compare to" TYPE STRING DEFAULT "CTRL" (Name of the sample group of which you want to identify the differentially expressed genes of.)
# PARAMETER cluster: "Name of the cluster" TYPE STRING DEFAULT 3 (Name of the cluster of which you want to identify the differentially expressed genes of. By default, the clusters are named with numbers starting from 0.)
# PARAMETER OPTIONAL test.type: "Which test to use for detecting marker genes" TYPE [wilcox, DESeq2, bimod, roc, t, tobit, poisson, negbinom] DEFAULT wilcox (Seurat currently implements Wilcoxon rank sum test, bimod \(likelihood-ratio test for single cell gene expression\), roc \(standard AUC classifier\), Students t-test, Tobit-test, poisson, negbinom and DESeq2. The latter three should be used on UMI datasets only, and assume an underlying poisson or negative-binomial distribution. Note that DESeq2 can sometimes give as results NAs for some genes.)
# PARAMETER OPTIONAL only.positive: "Return only positive marker genes" TYPE [FALSE, TRUE] DEFAULT TRUE (Tool only returns positive markers as default. Change the parameter here if you want to also include the negative markers.)
# PARAMETER OPTIONAL logFC.de: "Fold change threshold for differentially expressed genes in log2 scale" TYPE DECIMAL FROM 0 TO 5 DEFAULT 0.25 (Genes with an average fold change smaller than this are not included in the analysis.)
# PARAMETER OPTIONAL pval.cutoff.de: "Adjusted p-value cutoff for differentially expressed genes" TYPE DECIMAL FROM 0 TO 1 DEFAULT 0.05 (Cutoff for the adjusted p-value of the DE genes: by default, adjusted p-values bigger than 0.05 are filtered out.)
# PARAMETER OPTIONAL minpct: "Limit testing for differentially expressed genes to genes which are expressed in at least this fraction of cells" TYPE DECIMAL DEFAULT 0.1 (Test only genes which are detected in at least this fraction of cells in either of two samples being compared in the cluster of question. Meant to speed up testing by leaving out genes that are very infrequently expressed.)
# PARAMETER OPTIONAL replace.NAs: "How to handle DESeq2 NA values?" TYPE [replace, remove] DEFAULT remove (DESeq2 can sometimes give as results NAs for p-values for some genes. You can choose to remove them or replace NAs with 0.)
# RUNTIME R-4.3.2-single-cell
# SLOTS 2
# TOOLS_BIN ""

# 2024-03-07 ML

# Based on: https://satijalab.org/seurat/articles/parsebio_sketch_integration

library(Seurat)
library(dplyr)
library(tidyr)
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

# do this later?
if (normalisation.method == "SCT") {
  # When SCTransform was used to normalise the data, do a prep step:
  data.combined <- PrepSCTFindMarkers(data.combined)
}


# Rename clusters, because if they are just numbers, Seurat adds a "g" in front of the name,
# which can be confusing and looks stupid in the plots.
# data.combined$seurat_clusters[1:5]
data.combined$seurat_clusters <- paste("cluster", data.combined$seurat_clusters, sep="")

# Organise the levels alphabetically (as for some reason, they are not by default)
data.combined$stim <- as.ordered(data.combined$stim)
data.combined$type <- as.ordered(data.combined$type)
data.combined$seurat_clusters <- as.ordered(data.combined$seurat_clusters)
data.combined$orig.ident <- as.ordered(data.combined$orig.ident)

# Set the identity as clusters:
sel.clust <- "seurat_clusters"
data.combined <- SetIdent(data.combined, value = sel.clust)
# table(data.combined@active.ident)

# Aggregate 
# bulk <- AggregateExpression(data.combined,  return.seurat = T, assays = "RNA",  group.by = c("seurat_clusters", "stim", "type"), verbose = F) 
# ei assays = "RNA", so that it works also for sct?
bulk <- AggregateExpression(data.combined,  return.seurat = T, group.by = c("seurat_clusters", "stim", "type"), verbose = F) 

# Subset based on cluster number
cluster.bulk <- subset(bulk, seurat_clusters == paste("cluster",cluster, sep="")) 

# data.combined$seurat_clusters[1:5]
# cluster.bulk$seurat_clusters[1:5]

# Organise the levels alphabetically (as before):
cluster.bulk$stim <- as.ordered(cluster.bulk$stim)
cluster.bulk$type <- as.ordered(cluster.bulk$type)
cluster.bulk$seurat_clusters <- as.ordered(cluster.bulk$seurat_clusters)
cluster.bulk$orig.ident <- as.ordered(cluster.bulk$orig.ident)


# "type" contains the sample group information from the setup tool
Idents(cluster.bulk) <- "type"


# Compute differential expression
if (normalisation.method == "SCT") {
  #  save(cluster.bulk, file = "data_combined_pseudobulk.Robj")
  pseudobulk_markers <- FindMarkers(cluster.bulk, assay = "SCT", ident.1 = samples1, ident.2 = samples2, slot = "counts", , test.use = test.type, verbose = F, logfc.threshold = logFC.de, min.pct = minpct, recorrect_umi = FALSE, only.pos = only.positive)
  # Currently, it is unsure how to (best) use SCTransformed data with DESeq2:
  # https://github.com/satijalab/seurat/issues/7659
  if (test.type == "DESeq2"){
    stop("CHIPSTER-NOTE: Apologies, but currently in Seurat, DESeq2 won't work with SCTransformed data. Try with another method.")
  }
} else {
  pseudobulk_markers <- FindMarkers(cluster.bulk, ident.1 = samples1, ident.2 = samples2, slot = "counts", test.use = test.type, verbose = F, min.pct = minpct, only.pos = only.positive, logfc.threshold = logFC.de)
}

# FindMarkers function gives NAs to some "outlier" genes when using DESeq.
# https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#pvaluesNA
# "If a row contains a sample with an extreme count outlier then the p value and adjusted p value will be set to NA."" 

# Sepate NAs to another table for user:
pseudobulk_markers_NAs <- pseudobulk_markers[is.na(pseudobulk_markers$p_val_adj), ]
if(length(pseudobulk_markers_NAs) != 0){
 write.table(pseudobulk_markers_NAs, file = "pseudobulk_markers_NAs.tsv", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
}

# Remove or eplace NAs with 0.0, so that it doesn't cause problems later on.
if(replace.NAs == "replace"){
  pseudobulk_markers$p_val <- replace_na(pseudobulk_markers$p_val, 0.0)
  pseudobulk_markers$p_val_adj <- replace_na(pseudobulk_markers$p_val_adj, 0.0)
}else if (replace.NAs == "remove"){
  pseudobulk_markers <- pseudobulk_markers[!is.na(pseudobulk_markers$p_val_adj), ]
}

# Filter based on adj-p-val :
pseudobulk_markers_filtered <- pseudobulk_markers[pseudobulk_markers$p_val_adj < pval.cutoff.de, ]


## Comparison name for the output file:
comparison.name <- paste(samples1, "vs", samples2, "in_cluster", cluster, sep = "_")
name.for.output.file <- paste("pseudobulk_de-list_", comparison.name, ".tsv", sep = "")


# Plot
pdf(file = "expressionPlots.pdf")

pseudobulk_markers$gene <- rownames(pseudobulk_markers)
ggplot(pseudobulk_markers, aes(avg_log2FC, -log10(p_val))) + geom_point(size = 0.5, alpha = 0.5) + theme_bw() +
  ylab("-log10(unadjusted p-value)") + geom_text_repel(aes(label = ifelse(p_val_adj < 0.01, gene, "")), 
  colour = "red", size = 3)  + ggtitle("Volcano plot of pseudobulk results (p_val_adj < 0.01 in red)")

# Choose top 6, use head in case of ties
pseudobulk_markers %>%
  # group_by(cluster) %>%
  top_n(-6, p_val) %>%
  head(6)-> top5_cell_selection


# VlnPlot(cell_selection, features = as.character(unique(top5_cell_selection$gene)), ncol = 3, group.by = "type", assay = "RNA", pt.size = 0.1)
plotlist <- VlnPlot(cluster.bulk, features = as.character(rownames(top5_cell_selection)), ncol = 3, group.by = "type", assay = "RNA", pt.size = 0.1, combine = FALSE)
p1 <- cowplot::plot_grid(plotlist = plotlist, ncol = 2, label_size = 10)
title <- ggdraw() + draw_label("top genes, pseudobulk", fontface = 'bold')
cowplot::plot_grid(title, p1, ncol = 1, rel_heights = c(0.1, 1))


# Plot these genes across all clusters, but split by “type”
plotlist2 <- VlnPlot(data.combined, features = as.character(rownames(top5_cell_selection)), ncol = 3, split.by = "type", assay = "RNA", pt.size = 0.1, combine = FALSE) 
p2 <- cowplot::plot_grid(plotlist = plotlist2, ncol = 2, label_size = 10)
title <- ggdraw() + draw_label("all cells top genes, clusters", fontface = 'bold')
cowplot::plot_grid(title, p2, ncol = 1, rel_heights = c(0.1, 1))

Idents(data.combined) <- "stim"
plotlist3 <- VlnPlot(data.combined, features = as.character(rownames(top5_cell_selection)), ncol = 3, split.by = "type", assay = "RNA", pt.size = 0.1, combine = FALSE)  
p3 <- cowplot::plot_grid(plotlist = plotlist3, ncol = 2, label_size = 10)
title <- ggdraw() + draw_label("all cells top genes, samples", fontface = 'bold')
cowplot::plot_grid(title, p3, ncol = 1, rel_heights = c(0.1, 1))

# According: https://nbisweden.github.io/workshop-scRNAseq/labs/seurat/seurat_05_dge.html#pseudobulk
# Choose top 10, use head in case of ties
pseudobulk_markers %>%
  # group_by(cluster) %>%
  top_n(-10, p_val) %>%
  head(10) -> top10_cell_selection


#DotPlot(cluster.bulk,
 #       features = rownames(top10_cell_selection), group.by = "orig.ident",
#        assay = "RNA"
#) + coord_flip() + ggtitle("Expression, pseudobulk") + RotatedAxis()

# To get the samples in reasonable order:
dot_plot <- DotPlot(cluster.bulk, features = rownames(top10_cell_selection), group.by = "orig.ident", assay = "RNA")
levels_in_data <- unique(dot_plot$data$id)
dot_plot$data$id <- factor(x = dot_plot$data$id, levels = levels_in_data) # change the order of the factor levels
dot_plot + coord_flip() + ggtitle("Expression, pseudobulk") + RotatedAxis() # print again the plot

# For some reason, any of the following doesn't work:
# levels(as.ordered(as.character(dot_plot$data$id)))
# levels(gtools::mixedsort(dot_plot$data$id)))
# as.ordered(levels(dot_plot$data$id))
# order(as.character(levels(dot_plot$data$id)))
# sort(levels(dot_plot$data$id))
# order(levels(dot_plot$data$id))
# order(as.character(levels(dot_plot$data$id)))
# This works in R 4.2.3 (Seurat 5.0.1): levels(as.ordered(as.character(dot_plot$data$id)))

# save(dot_plot, file = "dot_plot.Robj")

DotPlot(data.combined,
        features = rownames(top10_cell_selection), group.by = "stim",
        assay = "RNA"
) + coord_flip() + ggtitle("Expression, all cells") + RotatedAxis()

dev.off() # close the pdf



 # Add average expression to the table:
   if (normalisation.method == "SCT") {
     # testaa viel tää: 
     aver_expr <- AverageExpression(object = cluster.bulk, slot = "data", assay = "SCT")
   } else {
     aver_expr <- AverageExpression(object = cluster.bulk)
   }

   aver_expr_in_clusters <- aver_expr[[1]]
   
   # select the wanted columns (based on samples1.cluster and samples2.cluster ) and rows (DEGs):
   aver_expr_ident1 <- round(aver_expr_in_clusters[row.names(pseudobulk_markers_filtered),samples1 ], digits = 4)
   aver_expr_ident2 <- round(aver_expr_in_clusters[row.names(pseudobulk_markers_filtered),samples2 ], digits = 4)

   full_table <- cbind(pseudobulk_markers_filtered, aver_expr_ident1 , aver_expr_ident2)
   colnames(full_table) <- c(colnames(full_table)[1:5], paste("aver_expr_",samples1), paste("aver_expr_",samples2))



# Comparison name for the output file:
 comparison.name <- paste(samples1, "vs", samples2, "in_cluster", cluster, sep = "_")
 name.for.output.file <- paste("pseudobulk_de-list_", comparison.name, ".tsv", sep = "")

# Write to table
 write.table(full_table, file = name.for.output.file, sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)



## EOF
