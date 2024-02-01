# TOOL single-cell-seurat-diffexp-chosen-samples.R: "Seurat v4 -Find DE genes between chosen samples" (This tool lists the differentially expressed genes between two user defined conditions or samples or sample groups for a user defined cluster. This tool can be used for Seurat objects with more than 2 samples in them.)
# INPUT OPTIONAL combined_seurat_obj.Robj: "Combined Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL de-list.tsv
# OUTPUT OPTIONAL de-list_{...}.tsv
# PARAMETER OPTIONAL normalisation.method: "Normalisation method used previously" TYPE [LogNormalize:"Global scaling normalization", SCT:"SCTransform"] DEFAULT LogNormalize (Which normalisation method was used in preprocessing, Global scaling normalization \(default, NormalizeData function used\) or SCTransform.)
# PARAMETER samples1: "Name of the samples to compare with" TYPE STRING DEFAULT "STIM" (Name of the sample or samples of which you want to identify the differentially expressed of.)
# PARAMETER samples2: "Name of the samples to compare to" TYPE STRING DEFAULT "CTRL" (Name of the sample or samples which you want to identify the differentially expressed of.)
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

# Load the R-Seurat-objects
load("combined_seurat_obj.Robj")

if (exists("seurat_obj")) {
  data.combined <- seurat_obj
}

DefaultAssay(data.combined) <- "RNA" # this is very crucial.

# Differentially expressed genes across conditions for the cluster (defined by the user, for example cluster 3 -> "3")
data.combined$celltype.stim <- paste(Idents(data.combined), data.combined$stim, sep = "_")
data.combined$celltype <- Idents(data.combined)
Idents(data.combined) <- "celltype.stim"

# Handle the sample listings:
samples1.ok <- unlist(strsplit(samples1, ", "))
samples2.ok <- unlist(strsplit(samples2, ", "))

# Check that the listed sample names are indeed in the data!
 if(sum(is.na(match(samples1.ok, data.combined$stim))) > 0) {
   which_one_missing <- samples1.ok[is.na(match(samples1.ok, data.combined$stim))]
   stop(print(paste("CHIPSTER-NOTE: Check the sample names given as parameter, sample name: '", which_one_missing, "' not found in the Seurat object.")))
 }

if(sum(is.na(match(samples2.ok, data.combined$stim))) > 0) {
  which_one_missing <- samples2.ok[is.na(match(samples2.ok, data.combined$stim))]
   stop(print(paste("CHIPSTER-NOTE: Check the sample names given as parameter, sample name: '", which_one_missing, "' not found in the Seurat object.")))
}


# Add the cluster name to the sample names:
samples1.cluster <- paste(cluster, "_", samples1.ok, sep = "")
samples2.cluster <- paste(cluster, "_", samples2.ok, sep = "")

# When SCTransform was used to normalise the data, do a prep step:
if (normalisation.method == "SCT") {
  data.combined <- PrepSCTFindMarkers(data.combined)
  # Note: assay = "SCT" and recorrect_umi = FALSE
  cluster_response <- FindMarkers(data.combined, assay = "SCT", ident.1 = samples1.cluster, ident.2 = samples2.cluster, verbose = FALSE, log2FC.threshold = logFC.de, min.pct = minpct, return.thresh = pval.cutoff.de, recorrect_umi = FALSE, only.pos = only.positive) 
} else {
  cluster_response <- FindMarkers(data.combined, ident.1 = samples1.cluster, ident.2 = samples2.cluster, verbose = FALSE, log2FC.threshold = logFC.de, min.pct = minpct, return.thresh = pval.cutoff.de, only.pos = only.positive)
}

 # Filter based on adj-p-val (no return.thresh parameter):
  cluster_response_filtered <- cluster_response[cluster_response$p_val_adj<pval.cutoff.de, ]

# Add average expression to the table:
if (normalisation.method == "SCT") {
  aver_expr <- AverageExpression(object = data.combined, slot = "data", assay = "SCT")
} else {
  aver_expr <- AverageExpression(object = data.combined)
}

aver_expr_in_clusters <- aver_expr[[1]]

# select the wanted columns (based on samples1.cluster and samples2.cluster ) and rows (DEGs):
aver_expr_ident1 <- round(aver_expr_in_clusters[row.names(cluster_response_filtered), samples1.cluster], digits = 4)
aver_expr_ident2 <- round(aver_expr_in_clusters[row.names(cluster_response_filtered), samples2.cluster], digits = 4)

full_table <- cbind(cluster_response_filtered, aver_expr_ident1, aver_expr_ident2)


# Comparison name for the output file:
comparison.name <- paste(samples1, "vs", samples2, "in_cluster", cluster, sep = "_")
name.for.output.file <- paste("de-list_", comparison.name, ".tsv", sep = "")

# Write to table
write.table(full_table, file = name.for.output.file, sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)

# Save the Robj for the next tool
# save(combined_seurat_obj, file="seurat_obj_combined.Robj")

## EOF
