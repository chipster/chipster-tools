# TOOL single-cell-seurat-diffexp-samples.R: "Seurat -Find conserved cluster markers and DE genes in two samples" (This tool lists the cell type markers that are conserved across the two conditions, and the differentially expressed genes between the two conditions for a user defined cluster. This tool can be used for two sample combined Seurat objects.) 
# INPUT OPTIONAL combined_seurat_obj.Robj: "Combined Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL de-list.tsv
# OUTPUT OPTIONAL conserved_markers.tsv
# PARAMETER cluster: "Name of the cluster" TYPE STRING DEFAULT 3 (Name of the cluster of which you want to identify the differentially expressed of. By default, the clusters are named with numbers starting from 0.)
# RUNTIME R-3.4.3



# 2018-16-05 ML

library(Seurat)

# Load the R-Seurat-objects (called seurat_obj -that's why we need to rename them here)
load("combined_seurat_obj.Robj")
#combined_seurat_obj <- data.combined

# Identify conserved cell type markers
# (uses package "metap" instead of metaDE since Seurat version 2.3.0)
nk.markers <- FindConservedMarkers(data.combined, ident.1 = cluster, grouping.var = "stim", 
		print.bar = FALSE)
#head(nk.markers)
write.table(nk.markers, file="conserved_markers.tsv", sep="\t", row.names=T, col.names=T, quote=F)


# Differentially expressed genes across conditions for the cluster (defined by the user, for example cluster 3 -> "3")
# cluster <- "3"
data.combined@meta.data$celltype.stim <- paste0(data.combined@ident, "_", 
		data.combined@meta.data$stim)
data.combined <- StashIdent(data.combined, save.name = "celltype")
data.combined <- SetAllIdent(data.combined, id = "celltype.stim")

lvls <- levels(as.factor(data.combined@meta.data$stim))
ident1 <- paste(cluster,"_", lvls[1], sep="")
ident2 <- paste(cluster,"_", lvls[2], sep="")
cluster_response <- FindMarkers(data.combined, ident.1 = ident1, ident.2 = ident2, 
		print.bar = FALSE)
# head(clustser_response, 15)
# show_rows <- 20
# cluster_response_rows <- head(cluster_response) #, show_rows)

write.table(cluster_response, file="de-list.tsv", sep="\t", row.names=T, col.names=T, quote=F)

# Save the Robj for the next tool
# save(combined_seurat_obj, file="seurat_obj_combined.Robj")

## EOF



