# TOOL deconv-test.R: "Deconv test" (Deconv test help text) 
# INPUT seurat_obj_subset.Robj: "Seurat object" TYPE GENERIC
# INPUT allen_cortex: "Allen Reference" TYPE GENERIC
# OUTPUT OPTIONAL seurat_obj_deconv.Robj
# OUTPUT OPTIONAL SpatialFeaturePlot.pdf
# OUTPUT OPTIONAL SpatialPlot.pdf
# OUTPUT OPTIONAL VlnPlot.pdf
# RUNTIME R-4.2.0-single-cell

library(dplyr, quietly = TRUE)
library(Seurat, quietly = TRUE)
suppressPackageStartupMessages(require(SCDC))
suppressPackageStartupMessages(require(Biobase))

system('ls')

# Load the reference dataset
allen_cortex <- readRDS("allen_cortex")

# Deconvolution:

allen_cortex@active.assay = "RNA"

markers_sc <- FindAllMarkers(allen_cortex, only.pos = TRUE, logfc.threshold = 0.1,
    test.use = "wilcox", min.pct = 0.05, min.diff.pct = 0.1, max.cells.per.ident = 200,
    return.thresh = 0.05, assay = "RNA")

markers_sc <- markers_sc[markers_sc$gene %in% rownames(seurat_obj),] # NEW 21.03

# Select top 20 genes per cluster, select top by first p-value, then absolute
# diff in pct, then quota of pct.
markers_sc$pct.diff <- markers_sc$pct.1 - markers_sc$pct.2
markers_sc$log.pct.diff <- log2((markers_sc$pct.1 * 99 + 1)/(markers_sc$pct.2 * 99 +
    1))
markers_sc %>%
    group_by(cluster) %>%
    top_n(-100, p_val) %>%
    top_n(50, pct.diff) %>%
    top_n(20, log.pct.diff) -> top20
m_feats <- unique(as.character(top20$gene))


eset_SC <- ExpressionSet(assayData = as.matrix(allen_cortex@assays$RNA@counts[m_feats,
    ]), phenoData = AnnotatedDataFrame(allen_cortex@meta.data))
eset_ST <- ExpressionSet(assayData = as.matrix(seurat_obj@assays$Spatial@counts[m_feats,
    ]), phenoData = AnnotatedDataFrame(seurat_obj@meta.data))

# Deconvolve

deconvolution_crc <- SCDC::SCDC_prop(bulk.eset = eset_ST, sc.eset = eset_SC, ct.varname = "subclass",
    ct.sub = as.character(unique(eset_SC$subclass)))


# head(deconvolution_crc$prop.est.mvw)

seurat_obj@assays[["SCDC"]] <- CreateAssayObject(data = t(deconvolution_crc$prop.est.mvw))

# Seems to be a bug in SeuratData package that the key is not set and any
# plotting function etc. will throw an error.
if (length(seurat_obj@assays$SCDC@key) == 0) {
    seurat_obj@assays$SCDC@key = "scdc_"
}

# Plotting:

pdf(file="SpatialFeaturePlot.pdf")

DefaultAssay(seurat_obj) <- "SCDC"
SpatialFeaturePlot(seurat_obj, features = c("L2/3 IT", "L4"), pt.size.factor = 1.6, ncol = 2,
    crop = TRUE)

dev.off()

pdf(file="SpatialPlot.pdf")

seurat_obj <- FindSpatiallyVariableFeatures(seurat_obj, assay = "SCDC", selection.method = "markvariogram",
    features = rownames(seurat_obj), r.metric = 5, slot = "data")
top.clusters <- head(SpatiallyVariableFeatures(seurat_obj), 4)
SpatialPlot(object = seurat_obj, features = top.clusters, ncol = 2)

dev.off()

pdf(file="VlnPlot.pdf")
VlnPlot(seurat_obj, group.by = "seurat_clusters", features = top.clusters, pt.size = 0,
    ncol = 2)
dev.off()

# 3.23:
# do I need dev.off() each time?
# get rid of all print statements
# find out why downloading of devtools is not working (prereq for SeuratData package)

# For the manual:
# Make it clear that:
# Add information about creating a NEW seurat subsetted object manually (not the one from the sample session) based
# on the clusters found in the UMAP 3rd figure (which was in the sample session) (pick 5 clusters using this 3rd figure)