suppressMessages(require(Seurat))
suppressMessages(require(scater))
suppressMessages(require(Matrix))

system("curl -O http://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v2/pbmc_1k_v2_filtered_feature_bc_matrix.h5")
system("curl -O http://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_filtered_feature_bc_matrix.h5")
system("curl -O http://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5")
v3.1k <- Read10X_h5("pbmc_1k_v3_filtered_feature_bc_matrix.h5", use.names = T)
v2.1k <- Read10X_h5("pbmc_1k_v2_filtered_feature_bc_matrix.h5", use.names = T)
p3.1k <- Read10X_h5("pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5", use.names = T)
# select only gene expression data
p3.1k <- p3.1k$`Gene Expression`

#create seurat objects
sdata.v2.1k <- CreateSeuratObject(v2.1k, project = "v2.1k")
sdata.v3.1k <- CreateSeuratObject(v3.1k, project = "v3.1k")
sdata.p3.1k <- CreateSeuratObject(p3.1k, project = "p3.1k")
# merge into one single seurat object. Add cell ids just in case you have overlapping barcodes between the
alldata <- merge(sdata.v2.1k, c(sdata.v3.1k,sdata.p3.1k), add.cell.ids=c("v2.1k","v3.1k","p3.1k"))
# also add in a metadata column that indicates v2 vs v3 chemistry
chemistry <- rep("v3",ncol(alldata))
chemistry[Idents(alldata) == "v2.1k"] <- "v2"
alldata <- AddMetaData(alldata, chemistry, col.name = "Chemistry")
alldata

# check number of cells from each sample, is stored in the orig.ident slot of metadata and is autmatically set as active ident.
table(Idents(alldata))

# Calculate mitochondrial proportion
mt.genes <- rownames(alldata)[grep("^MT-",rownames(alldata))]
C<-GetAssayData(object = alldata, slot = "counts")
percent.mito <- colSums(C[mt.genes,])/colSums(C)*100
alldata <- AddMetaData(alldata, percent.mito, col.name = "percent.mito")

# Calculate ribosomal proportion
rb.genes <- rownames(alldata)[grep("^RP[SL]",rownames(alldata))]
percent.ribo <- colSums(C[rb.genes,])/colSums(C)*100
alldata <- AddMetaData(alldata, percent.ribo, col.name = "percent.ribo")

# Plot QC
VlnPlot(alldata, features = "nFeature_RNA", pt.size = 0.1) + NoLegend()
VlnPlot(alldata, features = "nCount_RNA", pt.size = 0.1) + NoLegend()
VlnPlot(alldata, features = "percent.mito", pt.size = 0.1) + NoLegend()
VlnPlot(alldata, features = "percent.ribo", pt.size = 0.1) + NoLegend()
FeatureScatter(alldata, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(alldata, feature1 = "nFeature_RNA", feature2 = "percent.mito")
FeatureScatter(alldata, feature1="percent.ribo", feature2="nFeature_RNA")
FeatureScatter(alldata, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", cells = WhichCells(alldata, expression = orig.ident == "v3.1k") )
FeatureScatter(alldata, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", cells = WhichCells(alldata, expression = orig.ident == "v3.1k") )

# Filtering
#select cells with percent.mito < 25
selected <- WhichCells(alldata, expression = percent.mito < 25)
length(selected)
# and subset the object to only keep those cells
data.filt <- subset(alldata, cells = selected)
# plot violins for new data
VlnPlot(data.filt, features = 'percent.mito')

# gene detection filtering
#start with cells with many genes detected.
high.det.v3 <- WhichCells(data.filt, expression = nFeature_RNA > 4100)
high.det.v2 <- WhichCells(data.filt, expression = nFeature_RNA > 2000 & orig.ident == "v2.1k")
# remove these cells
data.filt <- subset(data.filt, cells=setdiff(WhichCells(data.filt),c(high.det.v2,high.det.v3)))
# check number of cells
ncol(data.filt)
# cells with low genes detected.
low.det.v3 <- WhichCells(data.filt, expression = nFeature_RNA < 1000 & orig.ident != "v2.1k")
low.det.v2 <- WhichCells(data.filt, expression = nFeature_RNA < 500 & orig.ident == "v2.1k")
# remove these cells
data.filt <- subset(data.filt, cells=setdiff(WhichCells(data.filt),c(low.det.v2,low.det.v3)))
# check number of cells
ncol(data.filt)

# Plot QC-stats again
VlnPlot(data.filt, features = "nFeature_RNA", pt.size = 0.1) + NoLegend()
VlnPlot(data.filt, features = "nCount_RNA", pt.size = 0.1) + NoLegend()
VlnPlot(data.filt, features = "percent.mito", pt.size = 0.1) + NoLegend()
VlnPlot(data.filt, features = "percent.ribo", pt.size = 0.1) + NoLegend()

# and check the number of cells per sample before and after filtering
table(Idents(alldata))
table(Idents(data.filt))

# Calculate cell-cycle scores
data.filt <- CellCycleScoring(
  object = data.filt,
  g2m.features = cc.genes$g2m.genes,
  s.features = cc.genes$s.genes
)
VlnPlot(data.filt, features = c("S.Score", "G2M.Score"))

# scater
sce <- as.SingleCellExperiment(data.filt)
# calculate all qc-metrics
sce <- calculateQCMetrics(sce, feature_controls = list(mito = mt.genes))
# check what all entries are - 
colnames(colData(sce))
colnames(rowData(sce))

## Plot QC stats

#Most expressed features
plotHighestExprs(sce, exprs_values = "counts")
# plot each sample separately
plotScater(sce, block1 = "ident", nfeatures = 1000)
# plot gene stats
plotRowData(sce, x = "n_cells_by_counts", y = "mean_counts")
#plot cell stats
p1 <- plotColData(sce, x = "total_counts", 
                  y = "total_features_by_counts", colour_by = "ident")
p2 <- plotColData(sce, x = "pct_counts_feature_control",
                  y = "total_features_by_counts", colour_by = "ident")
p3 <- plotColData(sce, x = "pct_counts_feature_control",
                  y = "pct_counts_in_top_50_features", colour_by = "ident")
multiplot(p1, p2, p3, cols = 2)

#identify outliers in qc-stats
sce <- runPCA(sce, use_coldata = TRUE, detect_outliers = TRUE)
plotReducedDim(sce, use_dimred="PCA_coldata", colour_by = "ident")
# check if we have any outliers
table(colData(sce)$outlier)

# run PCA with 1000 top variable genes
sce <- runPCA(sce, ntop = 1000, exprs_values = "logcounts", ncomponents = 20)
# PCA - with different coloring, first 4 components
# first by sample
plotPCA(sce,ncomponents=4,colour_by="ident")
# then by Celltype
plotPCA(sce,ncomponents=4,colour_by="percent.mito")

# Diffusion map, OBS! Requires installation of package destiny to run!
library(destiny)
set.seed(1)
sce <- runDiffusionMap(sce, ntop = 1000, ncomponents = 4)
plotDiffusionMap(sce, colour_by="ident",ncomponents=4)

# tSNE - uses Rtsne function to run tsne, here run with first 10 PCs
set.seed(1)
sce <- runTSNE(sce, ntop = 1000, ncomponents = 2, perplexity = 30, n_dimred = 10)
plotTSNE(sce, colour_by="ident")

# UMAP, OBS! Requires installation of package umap to run!
#library(umap)
set.seed(1)
sce <- runUMAP(sce)
plotUMAP(object = sce, colour_by="ident")

# explanatory factors
plotExplanatoryVariables(sce, variables =  c("ident","Chemistry","pct_counts_mito", "total_features_by_counts", "pct_counts_in_top_50_features", "pct_counts_in_top_500_features", "total_counts", "S.Score","G2M.Score"))
# for total features
plotExplanatoryPCs(sce, variables = c("ident", "Chemistry","pct_counts_mito", "total_features_by_counts", "pct_counts_in_top_50_features", "total_counts","S.Score","G2M.Score"), npcs_to_plot = 20)
sessionInfo()
