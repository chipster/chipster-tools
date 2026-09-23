# TOOL single-cell-seurat-de-between-samples-v5.R: "Seurat v5 -DE genes between samples for a given cluster" (For a given cluster, this tool finds differentially expressed genes by comparing each sample against all other samples. This tool can be used for combined Seurat objects with 2 or more samples.)
# INPUT OPTIONAL combined_seurat_obj.Robj: "Combined Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL de-list_{...}.tsv
# OUTPUT OPTIONAL log.txt
# PARAMETER OPTIONAL normalisation.method: "Normalisation method used previously" TYPE [LogNormalize:"Global scaling normalization", SCT:"SCTransform"] DEFAULT LogNormalize (Which normalisation method was used in preprocessing, Global scaling normalization \(default, NormalizeData function used\) or SCTransform.)
# PARAMETER cluster: "Name of the cluster" TYPE STRING DEFAULT 3 (Name of the cluster for which you want to find differentially expressed genes. By default, clusters are named with numbers starting from 0.)
# PARAMETER OPTIONAL only.positive: "Return only positive marker genes" TYPE [FALSE, TRUE] DEFAULT TRUE (Return only positive markers. Set to FALSE to also include negative markers.)
# PARAMETER OPTIONAL logFC.de: "Fold change in log2 scale" TYPE DECIMAL FROM 0 TO 5 DEFAULT 0.25 (Genes with an average fold change smaller than this are not included in the analysis.)
# PARAMETER OPTIONAL pval.cutoff.de: "Adjusted p-value cutoff" TYPE DECIMAL FROM 0 TO 1 DEFAULT 0.05 (Only genes with an adjusted p-value smaller than this are kept in the result table.)
# PARAMETER OPTIONAL minpct: "Limit testing to genes expressed in at least this fraction of cells" TYPE DECIMAL DEFAULT 0.1 (Test only genes detected in at least this fraction of cells in either of the two groups being compared. Speeds up testing by excluding very infrequently expressed genes.)
# PARAMETER OPTIONAL mincells: "Minimum number of cells in one of the groups" TYPE INTEGER DEFAULT 3 (If either the sample being tested or the combined group of all other samples has fewer than this many cells in the chosen cluster, that comparison is skipped and a warning is written to the log file.)
# RUNTIME R-4.3.2-single-cell
# SLOTS 5
# TOOLS_BIN ""


# 2026-09-21 ML Split off from single-cell-seurat-diffexp-samples-v5.R

library(Seurat)
options(Seurat.object.assay.version = "v5")

# Load the combined Seurat object
load("combined_seurat_obj.Robj")

if (exists("seurat_obj")) {
  data.combined <- seurat_obj
}

# When SCTransform was used, run the prep step first
if (normalisation.method == "SCT") {
  data.combined <- PrepSCTFindMarkers(data.combined)
}

DefaultAssay(data.combined) <- "RNA"

# Create per-cell labels combining cluster identity and sample name
data.combined$celltype.stim <- paste("cluster", Idents(data.combined), "-", data.combined$stim, sep = "")
data.combined$celltype <- Idents(data.combined)
Idents(data.combined) <- "celltype.stim"

lvls <- levels(as.factor(data.combined$stim))

# Sanity check: need at least 2 samples
if (length(lvls) < 2) {
  stop("CHIPSTER-NOTE: There are fewer than 2 samples in the data.")
}

# Open log file
log.file <- file("log.txt", open = "wt")
sink(log.file, append = TRUE, type = "output")
sink(log.file, append = TRUE, type = "message")

# For each sample: compare it against all others within the chosen cluster
for (i in 1:length(lvls)) {
  ident1 <- paste("cluster", cluster, "-", lvls[i], sep = "")
  ident2 <- paste("cluster", cluster, "-", lvls[-i], sep = "")

  # Check cell counts before running the test
  n_ident1 <- sum(Idents(data.combined) == ident1)
  n_ident2 <- sum(Idents(data.combined) %in% ident2)
  if (n_ident1 < mincells || n_ident2 < mincells) {
    message(paste0("Warning: skipping comparison '", lvls[i], " vs all others' for cluster ", cluster,
                   " — too few cells (", lvls[i], ": ", n_ident1, " cells; all others: ", n_ident2,
                   " cells). Minimum required: ", mincells, "."))
    next
  }

  if (normalisation.method == "SCT") {
    # Note: assay = "SCT" and recorrect_umi = FALSE required for SCTransform
    cluster_response <- FindMarkers(data.combined,
      assay = "SCT", ident.1 = ident1, ident.2 = ident2,
      verbose = FALSE, log2FC.threshold = logFC.de, min.pct = minpct,
      return.thresh = pval.cutoff.de, recorrect_umi = FALSE,
      only.pos = only.positive
    )
  } else {
    cluster_response <- FindMarkers(data.combined,
      ident.1 = ident1, ident.2 = ident2,
      verbose = FALSE, log2FC.threshold = logFC.de, min.pct = minpct,
      return.thresh = pval.cutoff.de, only.pos = only.positive
    )
  }

  # Filter on adjusted p-value
  cluster_response_filtered <- cluster_response[cluster_response$p_val_adj < pval.cutoff.de, ]

  # Add average expression columns for the compared groups
  # suppressMessages() is used here to silence Seurat v5's one-time nudge recommending  AggregateExpression() for pseudo-bulk analysis. That recommendation does not apply
  # here: AverageExpression() is used only to add informative average expression columns to the output table, not as a substitute for pseudo-bulk DE analysis.
  if (normalisation.method == "SCT") {
    aver_expr <- suppressMessages(AverageExpression(object = data.combined, slot = "data", assay = "SCT"))
} else {
    aver_expr <- suppressMessages(AverageExpression(object = data.combined))
}

  aver_expr_in_clusters <- aver_expr[[1]]
  aver_expr_ident1 <- round(aver_expr_in_clusters[row.names(cluster_response_filtered), ident1], digits = 4)
  aver_expr_ident2 <- round(aver_expr_in_clusters[row.names(cluster_response_filtered), ident2], digits = 4)

  full_table <- cbind(cluster_response_filtered, aver_expr_ident1, aver_expr_ident2)

  # Write one output file per sample
  comparison.name <- paste(lvls[i], "vsAllOthers", sep = "")
  name.for.output.file <- paste("de-list_", comparison.name, ".tsv", sep = "")
  write.table(full_table, file = name.for.output.file, sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
}

sink(type = "output")
sink(type = "message")
close(log.file)

# EOF