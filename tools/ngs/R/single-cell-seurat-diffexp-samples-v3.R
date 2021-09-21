# TOOL single-cell-seurat-diffexp-samples-v3.R: "Seurat v3 -Find conserved cluster markers and DE genes in two samples" (This tool lists the cell type markers that are conserved across the two conditions, and the differentially expressed genes between the two conditions for a user defined cluster. This tool can be used for two sample combined Seurat objects.) 
# INPUT OPTIONAL combined_seurat_obj.Robj: "Combined Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL de-list.tsv
# OUTPUT OPTIONAL conserved_markers.tsv
# PARAMETER cluster: "Name of the cluster" TYPE STRING DEFAULT 3 (Name of the cluster of which you want to identify the differentially expressed of. By default, the clusters are named with numbers starting from 0.)
# PARAMETER OPTIONAL only.positive: "Return only positive marker genes" TYPE [FALSE, TRUE] DEFAULT TRUE (Tool only returns positive markers as default. Change the parameter here if you want to also include the negative markers.)
# PARAMETER OPTIONAL logFC.conserved: "Fold change threshold for conserved markers in ln scale" TYPE DECIMAL FROM 0 TO 5 DEFAULT 0.25 (Genes with an average fold change smaller than this are not included in the analysis. Note that the base of log is e, so if you are interested in two-fold expression changes in linear scale, you need to enter 0.693 here.)
# PARAMETER OPTIONAL logFC.de: "Fold change threshold for differentially expressed genes in ln scale" TYPE DECIMAL FROM 0 TO 5 DEFAULT 0.25 (Genes with an average fold change smaller than this are not included in the analysis. Note that the base of log is e, so if you are interested in two-fold expression changes in linear scale, you need to enter 0.693 here.)
# PARAMETER OPTIONAL pval.cutoff.conserved: "Adjusted p-value cutoff for conserved markers" TYPE DECIMAL FROM 0 TO 1 DEFAULT 0.05 (Cutoff for the adjusted p-value of the conserved cluster marker genes: by default, adjusted p-values bigger than 0.05 are filtered out.)
# PARAMETER OPTIONAL pval.cutoff.de: "Adjusted p-value cutoff for differentially expressed genes" TYPE DECIMAL FROM 0 TO 1 DEFAULT 0.05 (Cutoff for the adjusted p-value of the DE genes: by default, adjusted p-values bigger than 0.05 are filtered out.)
# IMAGE comp-20.04-r-base
# RUNTIME R-4.1.0-single-cell


# RUNTIME R-3.6.1-single-cell


# 2018-16-05 ML 
# 11.07.2019 ML Seurat v3
# 23.09.2019 EK Add only.pos = TRUE
# 30.10.2019 ML Add filtering parameters
# 02.12.2019 EK Change FC filter to use logfc.threshold prefiltering, add adjusted p-value filtering for DE genes


library(Seurat)

# Load the R-Seurat-objects
load("combined_seurat_obj.Robj")

# Identify conserved cell type markers
# (uses package "metap" instead of metaDE since Seurat version 2.3.0)
DefaultAssay(data.combined) <- "RNA" # this is very crucial.
cluster.markers <- FindConservedMarkers(data.combined, ident.1 = cluster, grouping.var = "stim", only.pos = only.positive,
    verbose = FALSE, logfc.threshold = logFC.conserved)

# Filter based on logFC:
# PARAMETER logFC.cutoff.conserved: "Threshold for logFC of conserved markers" TYPE INTEGER DEFAULT 1 (Threshold for the logFC of the conserved cluster markers: by default, fold changes smaller than 1 are filtered out.)
# In case of negative fold changes;
# logFC.cutoff.conserved_2 <- -logFC.cutoff.conserved
# Note: hardcoded column names need to be changed to the stim levels in the next line
# dat2 <- subset(cluster.markers, (CTRL_avg_logFC >= logFC.cutoff.conserved | CTRL_avg_logFC <= logFC.cutoff.conserved_2) & (STIM_avg_logFC >= logFC.cutoff.conserved | STIM_avg_logFC <= logFC.cutoff.conserved_2)) 

# Filter conserved marker genes based on adj p-val:
# Note: hardcoded column names need to be changed to the stim levels
# dat2 <- subset(cluster.markers, (CTRL_p_val_adj < pval.cutoff.conserved & STIM_p_val_adj < pval.cutoff.conserved))
dat2 <- subset(cluster.markers, (cluster.markers[,5] < pval.cutoff.conserved & cluster.markers[,10] < pval.cutoff.conserved))

# Write to table
write.table(dat2, file = "conserved_markers.tsv", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)


# Differentially expressed genes across conditions for the cluster (defined by the user, for example cluster 3 -> "3")
data.combined$celltype.stim <- paste(Idents(data.combined), data.combined$stim, sep = "_")
data.combined$celltype <- Idents(data.combined)
Idents(data.combined) <- "celltype.stim"

lvls <- levels(as.factor(data.combined$stim))
ident1 <- paste(cluster, "_", lvls[1], sep = "")
ident2 <- paste(cluster, "_", lvls[2], sep = "")
cluster_response <- FindMarkers(data.combined, ident.1 = ident1, ident.2 = ident2, verbose = FALSE, logfc.threshold = logFC.de)

# Filter based on logFC:
# PARAMETER logFC.cutoff.de: "Threshold for logFC of DE genes" TYPE INTEGER DEFAULT 1 (Threshold for the logFC of the DE genes: by default, fold changes smaller than 1 are filtered out.)
# In case of negative fold changes;
# logFC.cutoff.de_2 <- -logFC.cutoff.de
# de2 <- subset(cluster_response, (avg_logFC >= logFC.cutoff.de | avg_logFC <= logFC.cutoff.de_2))

# Filter DE genes based on adj p-val:
de2 <- subset(cluster_response, (p_val_adj < pval.cutoff.de))

# Write to table
write.table(de2, file = "de-list.tsv", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)

# Save the Robj for the next tool
# save(combined_seurat_obj, file="seurat_obj_combined.Robj")

## EOF



