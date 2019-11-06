# TOOL single-cell-seurat-diffexp-samples-v3.R: "Seurat v3 -Find conserved cluster markers and DE genes in two samples" (This tool lists the cell type markers that are conserved across the two conditions, and the differentially expressed genes between the two conditions for a user defined cluster. This tool can be used for two sample combined Seurat objects.) 
# INPUT OPTIONAL combined_seurat_obj.Robj: "Combined Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL de-list.tsv
# OUTPUT OPTIONAL conserved_markers.tsv
# PARAMETER cluster: "Name of the cluster" TYPE STRING DEFAULT 3 (Name of the cluster of which you want to identify the differentially expressed of. By default, the clusters are named with numbers starting from 0.)
# PARAMETER OPTIONAL only.positive: "Only return positive markers" TYPE [FALSE, TRUE] DEFAULT TRUE (Tool only returns positive markers as default. Change the parameter here if you want to also include the negative markers.)
# PARAMETER logFC.cutoff.conserved: "Threshold for logFC of conserved markers" TYPE INTEGER DEFAULT 1 (Threshold for the logFC of the conserved cluster markers: by default, fold changes smaller than 1 are filtered out.)
# PARAMETER pval.cutoff.conserved: "Threshold for adjusted p-value of conserved markers" TYPE DECIMAL FROM 0 TO 1 DEFAULT 0.05 (Threshold for the adjusted p-value of the conserved cluster markers: by default, adjusted p-values bigger than 0.05 are filtered out.)
# RUNTIME R-3.6.1


# 2018-16-05 ML
# 11.07.2019 ML Seurat v3
# 2019-09-23 EK Add only.pos = TRUE
# 30.10.2019 ML Add filtering parameters

library(Seurat)

# Load the R-Seurat-objects
load("combined_seurat_obj.Robj")

# Identify conserved cell type markers
# (uses package "metap" instead of metaDE since Seurat version 2.3.0)
DefaultAssay(data.combined) <- "RNA" # this is very crucial.
nk.markers <- FindConservedMarkers(data.combined, ident.1 = cluster, grouping.var = "stim", only.pos = only.positive,
		verbose = FALSE)
#head(nk.markers)
# Filter based on logFC:
	# In case of negative fold changes;
	logFC.cutoff.conserved_2 <- -logFC.cutoff.conserved
  dat2 <- subset(nk.markers,(CTRL_avg_logFC>=logFC.cutoff.conserved | CTRL_avg_logFC<=logFC.cutoff.conserved_2) & (STIM_avg_logFC>=logFC.cutoff.conserved | STIM_avg_logFC<=logFC.cutoff.conserved_2))
# Filter based on adj p-val:
dat3 <- subset(dat2, (CTRL_p_val_adj<pval.cutoff.conserved & STIM_p_val_adj<pval.cutoff.conserved))


# Write to table:
write.table(dat3, file="conserved_markers.tsv", sep="\t", row.names=T, col.names=T, quote=F)


# Differentially expressed genes across conditions for the cluster (defined by the user, for example cluster 3 -> "3")
data.combined$celltype.stim <- paste(Idents(data.combined), data.combined$stim, sep = "_")
data.combined$celltype <- Idents(data.combined)
Idents(data.combined) <- "celltype.stim"

lvls <- levels(as.factor(data.combined$stim))
ident1 <- paste(cluster,"_", lvls[1], sep="")
ident2 <- paste(cluster,"_", lvls[2], sep="")
cluster_response <- FindMarkers(data.combined, ident.1 = ident1, ident.2 = ident2, verbose = FALSE)
# head(cluster.response, n = 15)

write.table(cluster_response, file="de-list.tsv", sep="\t", row.names=T, col.names=T, quote=F)

# Save the Robj for the next tool
# save(combined_seurat_obj, file="seurat_obj_combined.Robj")

## EOF



