# TOOL single-cell-seurat-clustering.R: "Seurat v2 -Clustering and detection of cluster marker genes" (Clusters cells, performs non-linear dimensional reduction tSNE for visualization purposes, and finds marker genes for the clusters.)
# INPUT seurat_obj.Robj: "Seurat object from the Setup tool" TYPE GENERIC
# OUTPUT OPTIONAL log.txt
# OUTPUT OPTIONAL tSNEplot.pdf
# OUTPUT OPTIONAL seurat_obj_2.Robj
# OUTPUT OPTIONAL markers.tsv
# PARAMETER OPTIONAL pcs_use: "Number of principal components to use" TYPE INTEGER DEFAULT 10 (How many principal components to use. User must define this based on the PCA-elbow and PCA plots from the setup tool.)
# PARAMETER OPTIONAL res: "Resolution for granularity" TYPE DECIMAL DEFAULT 0.6 (Resolution parameter that sets the granularity of the clustering. Increased values lead to greater number of clusters. Values between 0.6-1.2 return good results for single cell datasets of around 3K cells. For larger data sets, try higher resolution.)
# PARAMETER OPTIONAL perplex: "Perplexity, expected number of neighbors for tSNE plot" TYPE INTEGER DEFAULT 30 (Perplexity, expected number of neighbors. Default 30. Set to lower number if you have very few cells. Used for the tSNE visualisation of the clusters.)
# PARAMETER OPTIONAL point.size: "Point size in tSNE plot" TYPE DECIMAL DEFAULT 1 (Point size for tSNE plot.)
# PARAMETER OPTIONAL minpct: "Min fraction of cells where a cluster marker gene is expressed" TYPE DECIMAL DEFAULT 0.25 (Test only genes which are detected in at least this fraction of cells in either of the two populations. Meant to speed up the function by not testing genes that are very infrequently expression)
# PARAMETER OPTIONAL threshuse: "Differential expression threshold for a cluster marker gene" TYPE DECIMAL DEFAULT 0.25 (Limit testing to genes which show, on average, at least X-fold difference, in log2 scale, between the two groups of cells. Increasing thresh.use speeds up the function, but can miss weaker signals.)
# PARAMETER OPTIONAL test.type: "Which test to use for finding marker genes" TYPE [wilcox, bimod, roc, t, tobit, poisson, negbinom, MAST] DEFAULT wilcox (Denotes which test to use. Seurat currently implements \"wilcox\" \(Wilcoxon rank sum test, default\), \"bimod\" \(likelihood-ratio test for single cell gene expression\), \"roc\" \(standard AUC classifier\), \"t\" \(Students t-test\), \"tobit\" \(Tobit-test for differential gene expression\), \"MAST\" \(GLM-framework that treates cellular detection rate as a covariate\), \"poisson\", and \"negbinom\". The latter two options should be used on UMI datasets only, and assume an underlying poisson or negative-binomial distribution.)
# RUNTIME R-3.4.3



# 09.06.2017 ML
# 2018-01-11 ML update Seurat version to 2.2.0
# 2018-07-25 ML remove onlypos -parameter
# 2018-09-26 ML add perplexity parameter for datasets with fewer cells
# 2019-02-18 ML add heatmap plotting & number of cells in each cluster
# 2019-04-10 ML add more test.types (wilcox, MAST, DESeq2), wilcox changed as the default
# 2019-06-15 EK change name to include marker gene detection
# 2019-06-28 EK add point size parameter for tSNE plot in the code
# 2019-11-20 EK remove DESeq2 option

library(Seurat)
library(dplyr)
library(Matrix)
library(gplots)

point.size

# Load the R-Seurat-object (called seurat_obj)
load("seurat_obj.Robj")

#  Cluster the cells
# First set the directory for this step.
library(tools)
dir <- getwd()
dir <- file_path_as_absolute(dir)
# seurat_obj <- FindClusters(seurat_obj, pc.use=1:pcs_use, resolution = res, print.output= 0, save.SNN= T, temp.file.location = dir)
seurat_obj <- FindClusters(object = seurat_obj, reduction.type = "pca", dims.use = 1:pcs_use, 
		resolution = res, print.output = 0, save.SNN = TRUE)


# Non-linear dimensional reduction (tSNE) & number of cells in clusters
seurat_obj <- RunTSNE(seurat_obj, dims.use=1:pcs_use, do.fast=T, perplexity=perplex)

# Calculate number of cells per cluster from object@ident
cell.num <- table(seurat_obj@ident)

# Add cell number per cluster to cluster labels
ClusterLabels = paste("Cluster", names(cell.num), paste0("(n = ", cell.num, ")"))

# Order legend labels in plot in the same order as 'ClusterLabels'
ClusterBreaks = names(cell.num)

# Plot tSNE with new legend labels for clusters
pdf(file="tSNEplot.pdf") 
TSNEPlot(object = seurat_obj, do.return = T, pt.size = point.size, plot.title = paste("Number of cells: ", length(seurat_obj@cell.names))) +
		scale_colour_discrete(breaks = ClusterBreaks, 
				labels = ClusterLabels) +
		labs(x = "t-SNE 1",
				y = "t-SNE 2")

# Find all markers 
markers <- FindAllMarkers(seurat_obj, min.pct = minpct, thresh.use = threshuse, test.use = test.type) # min.pct = 0.25, thresh.use = 0.25, only.pos = onlypos

if(length(warnings())>0){ # or !is.null(warnings())
	stop("CHIPSTER-NOTE: There was issue with FindAllMarkers functions with the selected test type, try another test!")
}

write.table(as.matrix(markers), file = "markers.tsv", sep="\t", row.names=T, col.names=T, quote=F)

# Plot top10 genes of each cluster as a heatmap 
top10 <- markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
# setting slim.col.label to TRUE will print just the cluster IDS instead of
# every cell name
#### pdf(file="tSNEplot.pdf")
DoHeatmap(object = seurat_obj, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE, cex.row=4)

dev.off() # close the pdf


# Save the Robj for the next tool
save(seurat_obj, file="seurat_obj_2.Robj")

# EOF
