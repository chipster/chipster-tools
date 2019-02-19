# TOOL single-cell-seurat-clustering.R: "Seurat -Clustering" (Clusters cells, performs non-linear dimensional reduction tSNE for visualization purposes, and finds marker genes for the clusters.)
# INPUT seurat_obj.Robj: "Seurat object from the Setup tool" TYPE GENERIC
# OUTPUT OPTIONAL log.txt
# OUTPUT OPTIONAL tSNEplot.pdf
# OUTPUT OPTIONAL seurat_obj_2.Robj
# OUTPUT OPTIONAL markers.tsv
# PARAMETER OPTIONAL pcs_use: "Number of principal components to use" TYPE INTEGER DEFAULT 10 (How many principal components to use. User must define this based on the PCA-elbow and PCA plots from the setup tool.)
# PARAMETER OPTIONAL res: "Resolution for granularity" TYPE DECIMAL DEFAULT 0.6 (Resolution parameter that sets the granularity of the clustering. Increased values lead to greater number of clusters. Values between 0.6-1.2 return good results for single cell datasets of around 3K cells. For larger data sets, try higher resolution.)
# PARAMETER OPTIONAL perplex: "Perplexity, expected number of neighbors" TYPE INTEGER DEFAULT 30 (Perplexity, expected number of neighbors. Default 30. Set to lower number if you have very few cells.)
# PARAMETER OPTIONAL minpct: "Min fraction of cells where a cluster marker gene is expressed" TYPE DECIMAL DEFAULT 0.25 (Test only genes which are detected in at least this fraction of cells in either of the two populations. Meant to speed up the function by not testing genes that are very infrequently expression)
# PARAMETER OPTIONAL threshuse: "Differential expression threshold for a cluster marker gene" TYPE DECIMAL DEFAULT 0.25 (Limit testing to genes which show, on average, at least X-fold difference, in log-scale, between the two groups of cells. Increasing thresh.use speeds up the function, but can miss weaker signals.)
# PARAMETER OPTIONAL test.type: "Which test to use for finding marker genes" TYPE [bimod, roc, t, tobit, poisson, negbinom] DEFAULT bimod (Denotes which test to use. Seurat currently implements \"bimod\" \(likelihood-ratio test for single cell gene expression, McDavid et al., Bioinformatics, 2011, default\), \"roc\" \(standard AUC classifier\), \"t\" \(Students t-test\), and \"tobit\" \(Tobit-test for differential gene expression, as in Trapnell et al., Nature Biotech, 2014\), \"poisson\", and \"negbinom\". The latter two options should only be used on UMI datasets, and assume an underlying poisson or negative-binomial distribution.)
# RUNTIME R-3.4.3


# max.cells.per.ident ?
# min.pct = 0.25, thresh.use = 0.25
# removed: PARAMETER OPTIONAL onlypos: "Only positive changes" TYPE [TRUE, FALSE] DEFAULT FALSE (Only return positive markers if set to TRUE.)


# 09.06.2017 ML
# 2018-01-11 ML update Seurat version to 2.2.0
# 2018-07-25 ML removed onlypos -parameter
# 2018-09-26 ML add perplexity parameter for datasets with fewer cells
# 2019-02-18 ML add heatmap plotting

library(Seurat)
library(dplyr)
library(Matrix)
library(gplots)
library(tools)

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


# Non-linear dimensional reduction (tSNE)
seurat_obj <- RunTSNE(seurat_obj, dims.use=1:pcs_use, do.fast=T, perplexity=perplex)
pdf(file="tSNEplot.pdf") 
TSNEPlot(seurat_obj)

# Number of cells:
textplot(paste("\v \v Number of \n \v \v cells: \n \v \v", length(seurat_obj@cell.names)), halign="center", valign="center", cex=2) #, cex=0.8

# Find all markers 
markers <- FindAllMarkers(seurat_obj, min.pct = minpct, thresh.use = threshuse, test.use = test.type) # min.pct = 0.25, thresh.use = 0.25, only.pos = onlypos
write.table(as.matrix(markers), file = "markers.tsv", sep="\t", row.names=T, col.names=T, quote=F)

# Save the Robj for the next tool
save(seurat_obj, file="seurat_obj_2.Robj")

# Plot top10 genes of each cluster as a heatmap 
top10 <- markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
# setting slim.col.label to TRUE will print just the cluster IDS instead of
# every cell name
DoHeatmap(object = seurat_obj, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE)

dev.off() # close the pdf

# EOF
