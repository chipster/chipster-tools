# TOOL single-cell-seurat-diffexp-samples.R: "Seurat -Find conserved cluster markers and DE genes in two samples" (This tool lists the cell type markers that are conserved across the two conditions, and the differentially expressed genes between the two conditions for a user defined cluster. This tool can be used for two sample combined Seurat objects.) 
# INPUT OPTIONAL combined_seurat_obj.Robj: "Combined Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL de-list.tsv
# OUTPUT OPTIONAL conserved_markers.tsv
# PARAMETER cluster: "Name of the cluster" TYPE STRING DEFAULT 3 (Name of the cluster of which you want to identify the differentially expressed of. By default, the clusters are named with numbers starting from 0.)
# RUNTIME R-3.4.3


# 2018-16-05 ML
# 11.07.2019 ML Seurat v3

library(Seurat)

# Load the R-Seurat-objects
load("combined_seurat_obj.Robj")

# Identify conserved cell type markers
# (uses package "metap" instead of metaDE since Seurat version 2.3.0)
DefaultAssay(data.combined) <- "RNA"
nk.markers <- FindConservedMarkers(data.combined, ident.1 = cluster, grouping.var = "stim", 
		verbose = FALSE)
#head(nk.markers)
write.table(nk.markers, file="conserved_markers.tsv", sep="\t", row.names=T, col.names=T, quote=F)


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



