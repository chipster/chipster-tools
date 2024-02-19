# TOOL single-cell-seurat-diffexp-samples.R: "Seurat v4 -Find conserved cluster markers and DE genes in multiple samples" (This tool lists the cell type markers that are conserved across the conditions, and the differentially expressed genes between the conditions for a user defined cluster. In case of more than 2 samples, each sample is compared to all the other samples. This tool can be used for Seurat objects with 2 or more samples in them.)
# INPUT OPTIONAL combined_seurat_obj.Robj: "Combined Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL de-list.tsv
# OUTPUT OPTIONAL de-list_{...}.tsv
# OUTPUT OPTIONAL conserved_markers.tsv
# PARAMETER OPTIONAL normalisation.method: "Normalisation method used previously" TYPE [LogNormalize:"Global scaling normalization", SCT:"SCTransform"] DEFAULT LogNormalize (Which normalisation method was used in preprocessing, Global scaling normalization \(default, NormalizeData function used\) or SCTransform.)
# PARAMETER cluster: "Name of the cluster" TYPE STRING DEFAULT 3 (Name of the cluster of which you want to identify the differentially expressed of. By default, the clusters are named with numbers starting from 0.)
# PARAMETER OPTIONAL only.positive: "Return only positive marker genes" TYPE [FALSE, TRUE] DEFAULT TRUE (Tool only returns positive markers as default. Change the parameter here if you want to also include the negative markers.)
# PARAMETER OPTIONAL logFC.conserved: "Conserved markers: Fold change in log2 scale" TYPE DECIMAL FROM 0 TO 5 DEFAULT 0.25 (Genes with an average fold change smaller than this are not included in the analysis. )
# PARAMETER OPTIONAL pval.cutoff.conserved: "Conserved markers: p-value cutoff" TYPE DECIMAL FROM 0 TO 1 DEFAULT 0.05 (Cutoff for the p-value of the conserved cluster marker genes: by default, p-values bigger than 0.05 in any sample are filtered out.)
# PARAMETER OPTIONAL minpct_conserved: "Conserved markers: Limit testing to genes which are expressed in at least this fraction of cells" TYPE DECIMAL DEFAULT 0.1 (Test only genes which are detected in at least this fraction of cells in the cluster in question or in all the other cells. Meant to speed up testing by leaving out genes that are very infrequently expressed.)
# PARAMETER OPTIONAL mincellsconserved: "Conserved markers: Minimum number of cells in one of the groups" TYPE INTEGER DEFAULT 3 (How many cells at least there needs to be in each sample in the cluster in question.)
# PARAMETER OPTIONAL logFC.de: "Differentially expressed genes: Fold change in log2 scale" TYPE DECIMAL FROM 0 TO 5 DEFAULT 0.25 (Genes with an average fold change smaller than this are not included in the analysis.)
# PARAMETER OPTIONAL pval.cutoff.de: "Differentially expressed genes: Adjusted p-value cutoff" TYPE DECIMAL FROM 0 TO 1 DEFAULT 0.05 (Cutoff for the adjusted p-value of the DE genes: by default, adjusted p-values bigger than 0.05 are filtered out.)
# PARAMETER OPTIONAL minpct: "Differentially expressed genes: Limit testing to genes which are expressed in at least this fraction of cells" TYPE DECIMAL DEFAULT 0.1 (Test only genes which are detected in at least this fraction of cells in either of two samples being compared in the cluster of question. Meant to speed up testing by leaving out genes that are very infrequently expressed.)
# RUNTIME R-4.2.3-single-cell
# TOOLS_BIN ""



# # RUNTIME R-4.2.3-single-cell = old v4
# RUNTIME R-4.3.2-single-cell = new v5
# PARAMETER OPTIONAL returnthresh: "p-value threshold" TYPE DECIMAL DEFAULT 0.01 (Only return markers that have a p-value < return.thresh, or a power > return.thresh, if the test is ROC)


# 2018-16-05 ML
# 11.07.2019 ML Seurat v3
# 23.09.2019 EK Add only.pos = TRUE
# 30.10.2019 ML Add filtering parameters
# 02.12.2019 EK Change FC filter to use logfc.threshold prefiltering, add adjusted p-value filtering for DE genes
# 2021-10-04 ML Update to Seurat v4
# 2022-07-21 ML Tune for SCTransform data
# 2022-09-22 ML Fix conserved markers filtering for multiple sample case, add sanity check for cluster number
# 2022-09-20 ML Add min.cells.group parameter to allow outputting all the genes
# 2023-02-10 LG Add 2 slots
# 2023-04-03 ML Add parameters so that it is possible to print out all the genes + simplify code
# 2023-04-06 LG Remove 5 slots - Discrepancy in number of slots added v. removed


library(Seurat)
# options(Seurat.object.assay.version = "v5")

# Load the R-Seurat-objects
load("combined_seurat_obj.Robj")

if (exists("seurat_obj")) {
  data.combined <- seurat_obj
}

lvls <- levels(as.factor(data.combined$stim))

# When SCTransform was used to normalise the data, do a prep step:
if (normalisation.method == "SCT") {
  data.combined <- PrepSCTFindMarkers(data.combined)
}

# Check that this cluster is available in data:
# if is.na(match(cluster, levels(Idents(data.combined))) ){
#   stop("CHIPSTER-NOTE: Cluster given as input can not be found in the Seurat object!")
#  }

DefaultAssay(data.combined) <- "RNA" # this is very crucial.

# Identify conserved cell type markers
# (uses package "metap" instead of metaDE since Seurat version 2.3.0)
 if (normalisation.method == "SCT") {
    cluster.markers <- FindConservedMarkers(data.combined, assay = "SCT",
      ident.1 = cluster, grouping.var = "stim", only.pos = only.positive,
      verbose = FALSE, logfc.threshold = logFC.conserved, min.cells.group = mincellsconserved, min.pct = minpct_conserved, return.thresh = pval.cutoff.conserved)
  } else {
    cluster.markers <- FindConservedMarkers(data.combined,
      ident.1 = cluster, grouping.var = "stim", only.pos = only.positive,
      verbose = FALSE, logfc.threshold = logFC.conserved, min.cells.group = mincellsconserved, min.pct = minpct_conserved, return.thresh = pval.cutoff.conserved)
  }

# filter based on p_val_adj
# dat2 <- subset(cluster.markers, (CTRL_p_val_adj < pval.cutoff.conserved & STIM_p_val_adj < pval.cutoff.conserved))
# dat2 <- subset(cluster.markers, (cluster.markers[,5] < pval.cutoff.conserved & cluster.markers[,10] < pval.cutoff.conserved))
 dat2 <- subset(cluster.markers, cluster.markers[,"max_pval"] < pval.cutoff.conserved)

# install.packages("tidyverse", repos = "https://ftp.acc.umu.se/mirror/CRAN/") # remove this when the package is installed in tools!
# library("tidyverse")
# p.val.adj.table <- select(cluster.markers, ends_with("p_val_adj"))
# cluster.markers$max.adj.p.val <- apply(p.val.adj.table, 1, max, na.rm=TRUE)
# cluster.markers$minimum.adj.p.val <- apply(p.val.adj.table, 1, min, na.rm=TRUE)
# dat2 <- subset(cluster.markers, cluster.markers[,"max.adj.p.val"] < pval.cutoff.conserved)

# Write to table
write.table(dat2, file = "conserved_markers.tsv", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)


# Differentially expressed genes across conditions for the cluster (defined by the user, for example cluster 3 -> "3")
data.combined$celltype.stim <- paste("cluster", Idents(data.combined),"-", data.combined$stim, sep = "")
data.combined$celltype <- Idents(data.combined)
Idents(data.combined) <- "celltype.stim"


# Check that there are more than 2 samples:
if (length(lvls) < 2) {
  stop("CHIPSTER-NOTE: There are fewer than 2 samples in the data.")
} else {
  for (i in 1:length(lvls)) {
    # i:th sample vs all the others (= -i)
    ident1 <- paste("cluster", cluster, "-", lvls[i], sep = "")
    ident2 <- paste("cluster", cluster, "-", lvls[-i], sep = "")
    if (normalisation.method == "SCT") {
      # Note: assay = "SCT" and recorrect_umi = FALSE
      cluster_response <- FindMarkers(data.combined, assay = "SCT", ident.1 = ident1, ident.2 = ident2, verbose = FALSE, log2FC.threshold = logFC.de, min.pct = minpct, return.thresh = pval.cutoff.de, recorrect_umi = FALSE, only.pos = only.positive)
    } else {
      cluster_response <- FindMarkers(data.combined, ident.1 = ident1, ident.2 = ident2, verbose = FALSE, log2FC.threshold = logFC.de, min.pct = minpct, return.thresh = pval.cutoff.de, only.pos = only.positive)
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
   aver_expr_in_clusters[1:2,1:2]

   # select the wanted columns (based on samples1.cluster and samples2.cluster ) and rows (DEGs):
   aver_expr_ident1 <- round(aver_expr_in_clusters[row.names(cluster_response_filtered ), ident1], digits = 4)
   aver_expr_ident2 <- round(aver_expr_in_clusters[row.names(cluster_response_filtered ), ident2], digits = 4)

   full_table <- cbind(cluster_response_filtered , aver_expr_ident1, aver_expr_ident2)


   # Comparison name for the output file:
   comparison.name <- paste(lvls[i], "vsAllOthers", sep = "")
   name.for.output.file <- paste("de-list_", comparison.name, ".tsv", sep = "")

    # Write to table
    write.table(full_table, file = name.for.output.file, sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
  }
}
# Save the Robj for the next tool
# save(combined_seurat_obj, file="seurat_obj_combined.Robj")

## EOF
