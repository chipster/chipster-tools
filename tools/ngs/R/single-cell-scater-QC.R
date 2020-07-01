# TOOL single-cell-scater-QC.R: "Scater QC" (Quality control with Scater) 
# INPUT seurat_obj.Robj: "Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL log.txt
# OUTPUT OPTIONAL QCplots.pdf
# OUTPUT seurat_obj_2.Robj
# RUNTIME R-3.6.1


# 2020-06-23 ML

library(Seurat)
library(scater)

# Load the R-Seurat-object (called seurat_obj)
load("seurat_obj.Robj")

# pdf plots
pdf(file = "QCplots.pdf",, width = 13, height = 7)

sce <- as.SingleCellExperiment(seurat_obj)

mt.genes <- rownames(seurat_obj)[grep("^MT-",rownames(seurat_obj))]

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

# Not working atm, requires mvoutlier package that is not available for the current R version
#  https://github.com/ethering/scater_1.8.4_wrappers/issues/2 
# Identify outliers (=low quality cells): 
# run PCA and identify outliers in PCA space
# sce <- runPCA(sce, use_coldata = TRUE, detect_outliers = TRUE)
# plotReducedDim(sce, use_dimred="PCA_coldata") # colour_by = "ident"
# table(colData(sce)$outlier)
# textplot(paste("\v \v Number of \n \v \v outlier cells detected in PCA space: \n \v \v", table(colData(sce)$outlier)), halign = "center", valign = "center", cex = 2)

# run PCA with 1000 top variable genes
# Huom, ei tarvitsisi ajaa jos on jo Seurat-objectissa tämä, voisi testata enste.
sce <- runPCA(sce, ntop = 1000, exprs_values = "logcounts", ncomponents = 20)

# PCA - with different coloring, first 4 components
# color by sample
# plotPCA(sce,ncomponents=4,colour_by="ident")
# color by mito % and number of features:
plotPCA(sce,ncomponents=4,colour_by="percent.mt")
plotPCA(sce,ncomponents=4,colour_by="nFeature_RNA")
# 

# By default top 10 factors are plotted, but here we select some specific factors.
# Which factors to choose? Should we just use default?
# plotExplanatoryVariables(sce, variables =  c("ident","Chemistry","pct_counts_mito", "total_features_by_counts", "pct_counts_in_top_50_features", "pct_counts_in_top_500_features", "total_counts", "S.Score","G2M.Score"))
plotExplanatoryVariables(sce, variables =  c("pct_counts_mito", "total_features_by_counts", "pct_counts_in_top_50_features", "pct_counts_in_top_500_features", "total_counts"))
plotExplanatoryVariables(sce) 


# plotExplanatoryPCs(sce, variables = c("ident", "Chemistry","pct_counts_mito", "total_features_by_counts", "pct_counts_in_top_50_features", "total_counts","S.Score","G2M.Score"), npcs_to_plot = 20)
plotExplanatoryPCs(sce, variables =  c("pct_counts_mito", "total_features_by_counts", "pct_counts_in_top_50_features", "pct_counts_in_top_500_features", "total_counts"))
plotExplanatoryPCs(sce)

dev.off() # close the pdf

# Save the Robj for the next tool
save(seurat_obj, file = "seurat_obj_2.Robj")

# EOF