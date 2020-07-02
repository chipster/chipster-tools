# TOOL single-cell-scater-QC.R: "Scater QC" (Quality control with Scater. This tool takes as an input an Seurat object, and gives an pdf file with several quality control plots. The content of the plot depends on the contents of the Seurat object.) 
# INPUT seurat_obj.Robj: "Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL log.txt
# OUTPUT QCplots.pdf
# RUNTIME R-3.6.1

# 2020-06-23 ML

library(Seurat)
library(scater)
# library(mvoutlier)
library(gplots) # for textplot

# Load the R-Seurat-object (called seurat_obj)
load("seurat_obj.Robj")
# In case it's from the 2 sample analysis, change the name of the object:
# seurat_obj <- data.combined

# Transform Seurat object as Scater object:
sce <- as.SingleCellExperiment(seurat_obj)

# mito-genes:
mt.genes <- rownames(seurat_obj)[grep("^MT-",rownames(seurat_obj))]

# List available variables in the object:
available.variables <- colnames(colData(sce))

# Open pdf for plots
pdf(file = "QCplots.pdf",, width = 13, height = 7)

# Calculate all qc-metrics
sce <- calculateQCMetrics(sce, feature_controls = list(mito = mt.genes))

# Plot highest expression:
# Top 50 expressed genes. This can be valuable for detecting genes that are overabundant that may be driving a lot of the variation.
plotHighestExprs(sce, exprs_values = "counts", controls=NULL)

# Cumulative expression:
# Plot the relative proportion of the library size that is accounted for by the most highly expressed features for each cell (default 500 genes). 
# This can help us look for differences in expression distributions between samples.
plotScater(sce, nfeatures = 1000) # block1 = "ident"

# Plot cell stats:
p1 <- plotColData(sce, x = "total_counts", y = "total_features_by_counts") # colour_by = "ident"
p2 <- plotColData(sce, x = "pct_counts_feature_control", y = "total_features_by_counts")
p3 <- plotColData(sce, x = "pct_counts_feature_control",  y = "pct_counts_in_top_50_features")
p4 <- plotRowData(sce, x = "n_cells_by_counts", y = "mean_counts")
multiplot(p1, p2, p3, p4, cols = 2)

# Not recognizing library(mvoutlier) for some reason atm, so commented out:
# Identify outliers (=low quality cells): 
# Run PCA and identify outliers in PCA space
# sce <- runPCA(sce, use_coldata = TRUE, detect_outliers = TRUE)
# plotReducedDim(sce, use_dimred="PCA_coldata") # colour_by = "ident"
# table(colData(sce)$outlier)
# not_outliers <- table(colData(sce)$outlier)[1]
# outliers <- table(colData(sce)$outlier)[2]
# outliers_pt <- round(outliers/(outliers+not_outliers)*100, digits = 2)
# textplot(paste("Total number of cells", outliers+not_outliers, "\n Detecting outliers in PCA space: \n Not outlier cells:", not_outliers, " \n Outlier cells:", outliers, "(", outliers_pt ,"%)" ), halign = "center", valign = "center", cex = 1)
 
# Run PCA with 1000 top variable genes
# Note: not necessary, if Seurat object already has PCs? 
sce <- runPCA(sce, ntop = 1000, exprs_values = "logcounts", ncomponents = 20)

# PCA - with different coloring, first 4 components
# Color by sample:
# plotPCA(sce,ncomponents=4,colour_by="ident")
# color by mito % and number of features:
plotPCA(sce,ncomponents=4,colour_by="percent.mt")
plotPCA(sce,ncomponents=4,colour_by="nFeature_RNA")
# 

# Explanatory variables:
# By default top 10 factors are plotted, but we also select some specific factors.
# Top 10 factors:
plotExplanatoryVariables(sce) 
# If there are no cell-cycle scores available in the Seurat object:
if (is.na(pmatch("Phase",available.variables))) {
    plotExplanatoryVariables(sce, variables =  c("pct_counts_mito", "total_features_by_counts", "pct_counts_in_top_50_features", "pct_counts_in_top_500_features", "total_counts"))
    # plotExplanatoryVariables(sce, variables =  c("ident","Chemistry","pct_counts_mito", "total_features_by_counts", "pct_counts_in_top_50_features", "pct_counts_in_top_500_features", "total_counts", "S.Score","G2M.Score"))
}else {
   # If there are cell-cycle scores available, let's also plot them:
    plotExplanatoryVariables(sce, variables =  c("G2M.Score", "S.Score","Phase", "pct_counts_mito", "total_features_by_counts", "pct_counts_in_top_50_features", "pct_counts_in_top_500_features", "total_counts"))
}

# Explanatory PCs:
# Top 10 factors:
plotExplanatoryPCs(sce)
# If there are no cell-cycle scores available in the Seurat object:
if (is.na(pmatch("Phase",available.variables))) {
    # plotExplanatoryPCs(sce, variables = c("ident", "Chemistry","pct_counts_mito", "total_features_by_counts", "pct_counts_in_top_50_features", "total_counts","S.Score","G2M.Score"), npcs_to_plot = 20)
    plotExplanatoryPCs(sce, variables =  c("pct_counts_mito", "total_features_by_counts", "pct_counts_in_top_50_features", "pct_counts_in_top_500_features", "total_counts"))
}else {
   # If there are cell-cycle scores available, let's also plot them:
    plotExplanatoryPCs(sce, variables =  c("G2M.Score", "S.Score","Phase", "pct_counts_mito", "total_features_by_counts", "pct_counts_in_top_50_features", "pct_counts_in_top_500_features", "total_counts"))
}

dev.off() # close the pdf

# Not necessary: Save the Robj for the next tool
# OUTPUT seurat_obj_2.Robj
# save(seurat_obj, file = "seurat_obj_2.Robj")

# EOF