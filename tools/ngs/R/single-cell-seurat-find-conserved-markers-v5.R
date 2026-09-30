# TOOL single-cell-seurat-find-conserved-markers-v5.R: "Seurat v5 -Find conserved cluster markers in multiple samples" (For a given cluster, this tool finds cell type marker genes that are conserved across all conditions/samples. This tool can be used for combined Seurat objects with 2 or more samples.)
# INPUT OPTIONAL combined_seurat_obj.Robj: "Combined Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL conserved_markers.tsv
# OUTPUT OPTIONAL log.txt
# PARAMETER OPTIONAL normalisation.method: "Normalisation method used previously" TYPE [LogNormalize:"Global scaling normalization", SCT:"SCTransform"] DEFAULT LogNormalize (Which normalisation method was used in preprocessing, Global scaling normalization \(default, NormalizeData function used\) or SCTransform.)
# PARAMETER cluster: "Name of the cluster" TYPE STRING DEFAULT 3 (Name of the cluster for which you want to identify conserved marker genes. By default, clusters are named with numbers starting from 0.)
# PARAMETER OPTIONAL only.positive: "Return only positive marker genes" TYPE [FALSE, TRUE] DEFAULT TRUE (Return only positive markers. Set to FALSE to also include negative markers.)
# PARAMETER OPTIONAL logFC.conserved: "Fold change in log2 scale" TYPE DECIMAL FROM 0 TO 5 DEFAULT 0.25 (Genes with an average fold change smaller than this are not included in the analysis.)
# PARAMETER OPTIONAL pval.cutoff.conserved: "Adjusted p-value cutoff" TYPE DECIMAL FROM 0 TO 1 DEFAULT 0.05 (Only genes with an adjusted p-value smaller than this in all samples are kept in the result table.)
# PARAMETER OPTIONAL minpct_conserved: "Limit testing to genes expressed in at least this fraction of cells" TYPE DECIMAL DEFAULT 0.1 (Test only genes detected in at least this fraction of cells in the cluster in question or in all other cells. Speeds up testing by excluding very infrequently expressed genes.)
# PARAMETER OPTIONAL mincellsconserved: "Minimum number of cells in one of the groups" TYPE INTEGER DEFAULT 3 (Minimum number of cells required per sample in the cluster in question.)
# RUNTIME R-4.3.2-single-cell
# SLOTS 5
# TOOLS_BIN ""


# 2026-09-21 ML Split off from single-cell-seurat-diffexp-samples-v5.R

library(Seurat)
library(tidyverse)
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

# Identify conserved cell type markers across conditions
if (normalisation.method == "SCT") {
  cluster.markers <- FindConservedMarkers(data.combined,
    assay = "SCT",
    ident.1 = cluster, grouping.var = "stim", only.pos = only.positive,
    verbose = FALSE, logfc.threshold = logFC.conserved,
    min.cells.group = mincellsconserved, min.pct = minpct_conserved,
    return.thresh = pval.cutoff.conserved
  )
} else {
  cluster.markers <- FindConservedMarkers(data.combined,
    ident.1 = cluster, grouping.var = "stim", only.pos = only.positive,
    verbose = FALSE, logfc.threshold = logFC.conserved,
    min.cells.group = mincellsconserved, min.pct = minpct_conserved,
    return.thresh = pval.cutoff.conserved
  )
}

# Filter on maximum adjusted p-value across all samples
p.val.adj.table <- select(cluster.markers, ends_with("p_val_adj"))
cluster.markers$max.adj.p.val <- apply(p.val.adj.table, 1, max, na.rm = TRUE)
cluster.markers$minimum.adj.p.val <- apply(p.val.adj.table, 1, min, na.rm = TRUE)
dat2 <- subset(cluster.markers, cluster.markers[, "max.adj.p.val"] < pval.cutoff.conserved)

# Write to table
write.table(dat2, file = "conserved_markers.tsv", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)

# EOF