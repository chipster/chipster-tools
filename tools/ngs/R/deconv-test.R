# TOOL deconv-test.R: "Deconv test" (Deconv test help text) 

---
title: "Deconvolution"
output: html_document
date: "2023-02-13"
---

# Setup (with input files, libraries):

```{r fixing-Matrix-package-version-number-incompatibility-issue-with-R-version}
# packageVersion('Matrix') # 1.5.1
# suppressMessages(remove.packages('Matrix'))
# install.packages('Matrix', quiet = TRUE)
```

```{r fixing-Matrix-package-version-number-incompatibility-issue-with-R-version}
packageVersion('Matrix') # 1.5.3
library('Matrix')
# https://github.com/satijalab/seurat/issues/6746
```

```{r loading-libraries}
library(dplyr, quietly = TRUE)
library(Seurat, quietly = TRUE)
library(SeuratData, quietly = TRUE)
```

```{r loading-input-Seurat-object}
# load('/Users/lgerber/Downloads/seurat_obj_integrated.Robj')
# seurat_obj_integrated <- seurat_obj
# rm(seurat_obj)
```

```{r loading-input-Seurat-object}
# load('/Users/lgerber/Downloads/seurat_obj_subset.Robj')
# seurat_obj_subset <- seurat_obj
# rm(seurat_obj)
```

```{r loading-input-Seurat-object}
# load('/Users/lgerber/Downloads/seurat_spatial_obj_pca_INTEGRATED.Robj')
# seurat_spatial_obj_pca_INTEGRATED <- seurat_obj
# rm(seurat_obj)
```

```{r loading-input-Seurat-object}
load('/Users/lgerber/Downloads/seurat_spatial_obj_pca_ANTERIOR.Robj')
seurat_spatial_obj_pca_ANTERIOR <- seurat_obj
rm(seurat_obj)
```

```{r loading-input-allen-reference}
allen_reference <- readRDS("/Users/lgerber/Downloads/allen_cortex.rds")
```

# Subset ST for cortex_1

```{r}
# subset for the anterior dataset
cortex_1 <- subset(seurat_spatial_obj_pca_ANTERIOR, orig.ident == "anterior1")

# there seems to be an error in the subsetting, so the posterior1 image is not
# removed, do it manually
cortex_1@images$posterior1 = NULL

# subset for a specific region
# cortex_1 <- subset(cortex_1, anterior1_imagerow > 400 | anterior1_imagecol < 150, invert = TRUE)
# cortex_1 <- subset(cortex_1, anterior1_imagerow > 275 & anterior1_imagecol > 370, invert = TRUE)
# cortex_1 <- subset(cortex_1, anterior1_imagerow > 250 & anterior1_imagecol > 440, invert = TRUE)

# also subset for Frontal cortex_1 clusters
cortex_1 <- subset(cortex_1, idents = c(0, 1, 3, 6, 7))
# subset(cortex_1, idents = c(1, 2, 3, 4, 5)) # removed 5

p1 <- SpatialDimPlot(cortex_1, crop = TRUE, label = TRUE)
p2 <- SpatialDimPlot(cortex_1, crop = FALSE, label = TRUE, pt.size.factor = 1, label.size = 3)
p1 + p2
```


# Deconvolution:

```{r}
inst = installed.packages()

install.packages("remotes", quiet = TRUE)

if (!("xbioc" %in% rownames(inst))) {
    remotes::install_github("renozao/xbioc")
}
if (!("SCDC" %in% rownames(inst))) {
    remotes::install_github("meichendong/SCDC")
}

suppressPackageStartupMessages(require(SCDC))
suppressPackageStartupMessages(require(Biobase))
```

```{r}
allen_reference@active.assay = "RNA"

markers_sc <- FindAllMarkers(allen_reference, only.pos = TRUE, logfc.threshold = 0.1,
    test.use = "wilcox", min.pct = 0.05, min.diff.pct = 0.1, max.cells.per.ident = 200,
    return.thresh = 0.05, assay = "RNA")

markers_sc <- markers_sc[markers_sc$gene %in% rownames(cortex_1),] # NEW 21.03

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
```


```{r}
eset_SC <- ExpressionSet(assayData = as.matrix(allen_reference@assays$RNA@counts[m_feats,
    ]), phenoData = AnnotatedDataFrame(allen_reference@meta.data))
eset_ST <- ExpressionSet(assayData = as.matrix(cortex_1@assays$Spatial@counts[m_feats,
    ]), phenoData = AnnotatedDataFrame(cortex_1@meta.data))

```

# Deconvolve

```{r}
deconvolution_crc <- SCDC::SCDC_prop(bulk.eset = eset_ST, sc.eset = eset_SC, ct.varname = "subclass",
    ct.sub = as.character(unique(eset_SC$subclass)))
```

```{r}
head(deconvolution_crc$prop.est.mvw)
```

```{r}
cortex_1@assays[["SCDC"]] <- CreateAssayObject(data = t(deconvolution_crc$prop.est.mvw))

# Seems to be a bug in SeuratData package that the key is not set and any
# plotting function etc. will throw an error.
if (length(cortex_1@assays$SCDC@key) == 0) {
    cortex_1@assays$SCDC@key = "scdc_"
}
```

# Plotting:

```{r}
DefaultAssay(cortex_1) <- "SCDC"
SpatialFeaturePlot(cortex_1, features = c("L2/3 IT", "L4"), pt.size.factor = 1.6, ncol = 2,
    crop = TRUE)
# SpatialFeaturePlot(as.data.frame(cortex_1@assays[["SCDC"]]@data), features = c("L2/3 IT", "L4"), pt.size.factor = 1.6, ncol = 2, crop = TRUE)
                   
# as.data.frame(seurat_obj_integrated@assays[["SCDC"]]@data) # ncol = 2
```

```{r}
cortex_1 <- FindSpatiallyVariableFeatures(cortex_1, assay = "SCDC", selection.method = "markvariogram",
    features = rownames(cortex_1), r.metric = 5, slot = "data")
top.clusters <- head(SpatiallyVariableFeatures(cortex_1), 4)
SpatialPlot(object = cortex_1, features = top.clusters, ncol = 2)
```

```{r}
VlnPlot(cortex_1, group.by = "seurat_clusters", features = top.clusters, pt.size = 0,
    ncol = 2)
```
