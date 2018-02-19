# TOOL single-cell-seurat-clustering.R: "BETA Seurat -Clustering" (Clusters cells, performs non-linear dimensional reduction tSNE for visualization purposes, and finds marker genes for the clusters.)
# INPUT seurat_obj.Robj: "Seurat object from the Setup tool" TYPE GENERIC
# OUTPUT OPTIONAL log.txt
# OUTPUT OPTIONAL tSNEplot.pdf
# OUTPUT OPTIONAL seurat_obj_2.Robj
# OUTPUT OPTIONAL markers.tsv
# PARAMETER OPTIONAL pcs_use: "Number of principal components to use" TYPE INTEGER DEFAULT 10 (How many principal components to use. User must define this based on the PCA-elbow and PCA plots from the setup tool.)
# PARAMETER OPTIONAL res: "Resolution for granularity" TYPE DECIMAL DEFAULT 0.6 (Resolution parameter that sets the granularity of the clustering. Increased values lead to greater number of clusters. Values between 0.6-1.2 return good results for single cell datasets of around 3K cells. For larger data sets, try higher resolution.)
# PARAMETER OPTIONAL minpct: "Min fraction of cells where a cluster marker gene is expressed" TYPE DECIMAL DEFAULT 0.25 (Test only genes which are detected in at least this fraction of cells in either of the two populations. Meant to speed up the function by not testing genes that are very infrequently expression)
# PARAMETER OPTIONAL threshuse: "Differential expression threshold for a cluster marker gene" TYPE DECIMAL DEFAULT 0.25 (Limit testing to genes which show, on average, at least X-fold difference, in log-scale, between the two groups of cells. Increasing thresh.use speeds up the function, but can miss weaker signals.)
# PARAMETER OPTIONAL test.type: "Which test to use for finding marker genes" TYPE [bimod, roc, t, tobit, poisson, negbinom] DEFAULT bimod (Denotes which test to use. Seurat currently implements \"bimod\" \(likelihood-ratio test for single cell gene expression, McDavid et al., Bioinformatics, 2011, default\), \"roc\" \(standard AUC classifier\), \"t\" \(Students t-test\), and \"tobit\" \(Tobit-test for differential gene expression, as in Trapnell et al., Nature Biotech, 2014\), \"poisson\", and \"negbinom\". The latter two options should only be used on UMI datasets, and assume an underlying poisson or negative-binomial distribution.)
# PARAMETER OPTIONAL onlypos: "Only positive changes" TYPE [TRUE, FALSE] DEFAULT FALSE (Only return positive markers if set to TRUE.)
# RUNTIME R-3.4.3


# max.cells.per.ident ?
# min.pct = 0.25, thresh.use = 0.25


# 09.06.2017 ML
# 2018-01-11 ML update Seurat version to 2.2.0

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
seurat_obj <- RunTSNE(seurat_obj, dims.use=1:pcs_use, do.fast=T)
pdf(file="tSNEplot.pdf") 
TSNEPlot(seurat_obj)
dev.off() # close the pdf

# Find all markers 
markers <- FindAllMarkers(seurat_obj, only.pos = onlypos, min.pct = minpct, thresh.use = threshuse, test.use = test.type) # min.pct = 0.25, thresh.use = 0.25
write.table(as.matrix(markers), file = "markers.tsv", sep="\t", row.names=T, col.names=T, quote=F)

# Save the Robj for the next tool
save(seurat_obj, file="seurat_obj_2.Robj")


# EOF
