# TOOL single-cell-seurat-diffexp-chosen-samples.R: "Seurat v4 -Find DE genes between chosen samples" (This tool lists the differentially expressed genes between two user defined conditions or samples or sample groups for a user defined cluster. This tool can be used for Seurat objects with more than 2 samples in them.) 
# INPUT OPTIONAL combined_seurat_obj.Robj: "Combined Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL de-list.tsv
# OUTPUT OPTIONAL de-list_{...}.tsv
# PARAMETER samples1: "Name of the samples to compare to" TYPE STRING DEFAULT "CTRL" (Name of the sample or samples of which you want to identify the differentially expressed of.)
# PARAMETER samples2: "Name of the samples to compare with" TYPE STRING DEFAULT "STIM" (Name of the sample or samples which you want to identify the differentially expressed of.)
# PARAMETER cluster: "Name of the cluster" TYPE STRING DEFAULT 3 (Name of the cluster of which you want to identify the differentially expressed of. By default, the clusters are named with numbers starting from 0.)
# PARAMETER OPTIONAL only.positive: "Return only positive marker genes" TYPE [FALSE, TRUE] DEFAULT TRUE (Tool only returns positive markers as default. Change the parameter here if you want to also include the negative markers.)
# PARAMETER OPTIONAL logFC.de: "Fold change threshold for differentially expressed genes in log scale" TYPE DECIMAL FROM 0 TO 5 DEFAULT 0.25 (Genes with an average fold change smaller than this are not included in the analysis.)
# PARAMETER OPTIONAL pval.cutoff.de: "Adjusted p-value cutoff for differentially expressed genes" TYPE DECIMAL FROM 0 TO 1 DEFAULT 0.05 (Cutoff for the adjusted p-value of the DE genes: by default, adjusted p-values bigger than 0.05 are filtered out.)
# RUNTIME R-4.1.0-single-cell


# 2018-16-05 ML 
# 11.07.2019 ML Seurat v3
# 23.09.2019 EK Add only.pos = TRUE
# 30.10.2019 ML Add filtering parameters
# 02.12.2019 EK Change FC filter to use logfc.threshold prefiltering, add adjusted p-value filtering for DE genes
# 2021-10-04 ML Update to Seurat v4

library(Seurat)

# Load the R-Seurat-objects
load("combined_seurat_obj.Robj")

if (exists("seurat_obj")) {
  data.combined <- seurat_obj
}

# Differentially expressed genes across conditions for the cluster (defined by the user, for example cluster 3 -> "3")
data.combined$celltype.stim <- paste(Idents(data.combined), data.combined$stim, sep = "_")
data.combined$celltype <- Idents(data.combined)
Idents(data.combined) <- "celltype.stim"

# Handle the sample listings:
samples1.ok <- unlist(strsplit(samples1, ", "))
samples2.ok <- unlist(strsplit(samples2, ", "))
# Add the cluster name to the sample names:
samples1.cluster <- paste(cluster, "_", samples1.ok, sep="")
samples2.cluster <- paste(cluster, "_", samples2.ok, sep="")

cluster_response <- FindMarkers(data.combined, ident.1 = samples1.cluster, ident.2 = samples2.cluster, verbose = FALSE, logfc.threshold = logFC.de)
# Filter DE genes based on adj p-val:
de2 <- subset(cluster_response, (p_val_adj < pval.cutoff.de))

# Comparison name for the output file:
comparison.name <- paste(samples1, "vs", samples2, "in_cluster", cluster, sep="_")
name.for.output.file <- paste("de-list_", comparison.name, ".tsv", sep="")

# Write to table
write.table(de2, file = name.for.output.file, sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)

# Save the Robj for the next tool
# save(combined_seurat_obj, file="seurat_obj_combined.Robj")

## EOF



