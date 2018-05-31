# TOOL single-cell-seurat-integrated-analysis.R: "BETA Seurat -Integrated analysis" (This tool aligns the CCA subspaces and performs integrated analysis on the data.) 
# INPUT OPTIONAL combined_seurat_obj.Robj: "Combined Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL integrated_plot.pdf
# OUTPUT OPTIONAL seurat_obj_combined.Robj
# PARAMETER OPTIONAL num.dims: "Number of CCs to use " TYPE INTEGER DEFAULT 20 (Number of canonical correlates to use. Use the plots from CCA tool to estimate how many you wish to use.)
# PARAMETER OPTIONAL res: "Resolution for granularity" TYPE DECIMAL DEFAULT 0.6 (Resolution parameter that sets the granularity of the clustering. Increased values lead to greater number of clusters. Values between 0.6-1.2 return good results for single cell datasets of around 3K cells. For larger data sets, try higher resolution.)
# RUNTIME R-3.4.3




# PARAMETER OPTIONAL yhighcutoff: "Top cutoff on y-axis for identifying variable genes" TYPE DECIMAL DEFAULT Inf (For limiting the selection of variable genes.)
# OUTPUT OPTIONAL log.txt
# OUTPUT OPTIONAL de-list.tsv


# 2018-16-05 ML

library(Seurat)

# Load the R-Seurat-objects (called seurat_obj -that's why we need to rename them here)
load("combined_seurat_obj.Robj")
#combined_seurat_obj <- data.combined

# Align CCA subspaces
data.combined <- AlignSubspace(data.combined, reduction.type = "cca", grouping.var = "stim", 
		dims.align = 1:num.dims)

# We can visualize the aligned CCA and perform an integrated analysis.
pdf(file="integrated_plot.pdf", , width=13, height=7)  # open pdf
p1 <- VlnPlot(object = data.combined, features.plot = "ACC1", group.by = "stim", 
		do.return = TRUE)
p2 <- VlnPlot(object = data.combined, features.plot = "ACC2", group.by = "stim", 
		do.return = TRUE)
plot_grid(p1, p2)

# t-SNE and Clustering
data.combined <- RunTSNE(data.combined, reduction.use = "cca.aligned", dims.use = 1:num.dims, 
		do.fast = T)
data.combined <- FindClusters(data.combined, reduction.type = "cca.aligned", 
		resolution = res, dims.use = 1:num.dims)

# Visualization
p1 <- TSNEPlot(data.combined, do.return = T, pt.size = 0.5, group.by = "stim")
p2 <- TSNEPlot(data.combined, do.label = T, do.return = T, pt.size = 0.5)
plot_grid(p1, p2)


dev.off() # close the pdf


## Differentially expressed genes across conditions -pitää antaa mistä klusterista tämä halutaan.
#cluster <- "3"
#data.combined@meta.data$celltype.stim <- paste0(data.combined@ident, "_", 
#		data.combined@meta.data$stim)
#data.combined <- StashIdent(data.combined, save.name = "celltype")
#data.combined <- SetAllIdent(data.combined, id = "celltype.stim")
#
#lvls <- levels(as.factor(data.combined@meta.data$stim))
#ident1 <- paste(cluster,"_", lvls[1], sep="")
#ident2 <- paste(cluster,"_", lvls[2], sep="")
#cluster_response <- FindMarkers(data.combined, ident.1 = ident1, ident.2 = ident2, 
#		print.bar = FALSE)
## head(clustser_response, 15)
## show_rows <- 20
#cluster_response_rows <- head(cluster_response) #, show_rows)
#
## write.table(cbind(sig[,1:ndat], round(sig[, (ndat+1):(nmax-2)], digits=2), format(sig[, (nmax-1):nmax], digits=4, scientific=T)), file="de-list-deseq2.tsv", sep="\t", row.names=T, col.names=T, quote=F)
#write.table(cluster_response_rows, file="de-list.tsv", sep="\t", row.names=T, col.names=T, quote=F)
#
#diff_exp.tsv

# Save the Robj for the next tool
save(data.combined, file="seurat_obj_combined.Robj")

## EOF



