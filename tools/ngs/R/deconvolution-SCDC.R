# TOOL deconvolution-SCDC.R: "Seurat v4 and SCDC -Identify cell types using deconvolution" (This tool estimates the abundance \(or proportion\) of different cell types in Visium spots using a single cell reference dataset.)
# INPUT seurat_obj_subset.Robj: "Seurat Object Containing Visium Data" TYPE GENERIC
# INPUT allen_cortex: "Allen Reference" TYPE GENERIC
# OUTPUT OPTIONAL seurat_obj_deconv.Robj
# OUTPUT OPTIONAL SpatialFeaturePlot.pdf
# OUTPUT OPTIONAL SpatialPlot.pdf
# OUTPUT OPTIONAL VlnPlot.pdf
# PARAMETER top.genes: "Number of marker genes computed for each cell type in the reference dataset" TYPE INTEGER FROM 20 TO 100 DEFAULT 20 (Select the number \(between 20 and 100\) of marker genes \(for each cluster in the reference dataset\) to be used in the deconvolution step.)
# PARAMETER clusters: "Cell type(s) from the reference dataset to be plotted in SpatialFeaturePlot" TYPE STRING DEFAULT "L4" (If you list multiple cell types, please use comma\(s\) \(,\) as a separator, e.g., \"L2/3 IT\,L4\".)
# PARAMETER top_num_cluster: "Number of most spatially variable cell types to be plotted in SpatialPlot and violin plots" TYPE INTEGER FROM 1 TO 10 DEFAULT 4 (Select the number of top cell type\(s\) \(between 1 and 10\) to include in the analysis. These cell types are spatially restricted to particular location\(s\), as determined by the Markvariogram method.)
# RUNTIME R-4.2.0-single-cell
# SLOTS 4

library(dplyr, quietly = TRUE)
library(Seurat, quietly = TRUE)
# suppressPackageStartupMessages(require(SCDC))
suppressPackageStartupMessages(require(Biobase))
# library(githubinstall)

# install_github("meichendong/SCDC", ref = github_pull("31"))
packageDescription('SCDC')
system('ls')

# Load the reference dataset
allen_cortex <- readRDS("allen_cortex")

# Load the R-Seurat-object (called seurat_obj)
load("seurat_obj_subset.Robj")

clusters <- trimws(unlist(strsplit(clusters, ",")))

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
    top_n(top.genes, log.pct.diff) -> top_num_genes
m_feats <- unique(as.character(top_num_genes$gene))


eset_SC <- ExpressionSet(assayData = as.matrix(allen_cortex@assays$RNA@counts[m_feats,
    ]), phenoData = AnnotatedDataFrame(allen_cortex@meta.data))
eset_ST <- ExpressionSet(assayData = as.matrix(seurat_obj@assays$Spatial@counts[m_feats,
    ]), phenoData = AnnotatedDataFrame(seurat_obj@meta.data))

# Deconvolve

deconvolution_crc <- suppressMessages(SCDC::SCDC_prop(bulk.eset = eset_ST, sc.eset = eset_SC, ct.varname = "subclass",
    ct.sub = as.character(unique(eset_SC$subclass))))


# head(deconvolution_crc$prop.est.mvw)

seurat_obj@assays[["SCDC"]] <- CreateAssayObject(data = t(deconvolution_crc$prop.est.mvw))

# Seems to be a bug in SeuratData package that the key is not set and any
# plotting function etc. will throw an error.
if (length(seurat_obj@assays$SCDC@key) == 0) {
    seurat_obj@assays$SCDC@key = "scdc_"
}

# testing b/c Seurat object has incorrect assay
# seurat_obj@active.assay = 'SCDC'

save(seurat_obj, file="seurat_obj_deconv.Robj")

# Plotting:

pdf(file="SpatialFeaturePlot.pdf")

DefaultAssay(seurat_obj) <- "SCDC"
SpatialFeaturePlot(seurat_obj, features = clusters, pt.size.factor = 1.6, ncol = 2,
    crop = TRUE)

dev.off()


pdf(file="SpatialPlot.pdf")

seurat_obj <- FindSpatiallyVariableFeatures(seurat_obj, assay = "SCDC", selection.method = "markvariogram",
    features = rownames(seurat_obj), r.metric = 5, slot = "data")
top.clusters <- head(SpatiallyVariableFeatures(seurat_obj), top_num_cluster)
SpatialPlot(object = seurat_obj, features = top.clusters, ncol = 2)

dev.off()

pdf(file="VlnPlot.pdf")
VlnPlot(seurat_obj, group.by = "seurat_clusters", features = top.clusters, pt.size = 0,
    ncol = 2)
dev.off()