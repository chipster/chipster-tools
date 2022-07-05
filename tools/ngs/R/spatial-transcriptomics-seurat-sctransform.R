# TOOL spatial-transcriptomics-seurat-sctransform.R: "Seurat v4 -SCTransform: Filter cells, normalize, regress and detect high-variance features in spatial data" (This tool filters out dead cells, empties and doublets. It then normalizes gene expression values using the SCTransform method, detects highly variable genes, scales the data and regresses out unwanted variation based on the number of UMIs and mitochondrial transcript percentage. You can also choose to regress out variation due to cell cycle heterogeneity.)
# INPUT OPTIONAL seurat_spatial_setup.Robj: "Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL counts_plot.pdf 
# OUTPUT OPTIONAL seurat_obj_sctransform.Robj
# PARAMETER OPTIONAL mingenes: "Filter out cells which have less than this many genes expressed" TYPE INTEGER DEFAULT 500 (Filter out empties. The cells to be kept must express at least this number of genes.)
# PARAMETER OPTIONAL genecountcutoff: "Filter out cells which have more than this many genes expressed" TYPE INTEGER DEFAULT 2500 (Filter out multiplets. The cells to be kept must express less than this number of genes.)
# PARAMETER OPTIONAL mitocutoff: "Filter out cells which have higher mitochondrial transcript percentage" TYPE DECIMAL FROM 0 TO 100 DEFAULT 20 (Filter out dead cells. The cells to be kept must have lower percentage of mitochondrial transcripts than this.)
# PARAMETER OPTIONAL num.features: "Number of variable genes to return" TYPE INTEGER DEFAULT 3000 (Number of features to select as top variable features, i.e. how many features returned. For SCTransform, the recommended default is 3000.)
# RUNTIME R-4.1.0-single-cell

library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)

# Load the R-Seurat-object (called seurat_obj)
load("seurat_spatial_setup.Robj")

# Open the pdf file for plotting
pdf(file="counts_plot.pdf", width=13, height=7) 

plot1 <- VlnPlot(seurat_obj, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(seurat_obj, features = "nCount_Spatial") + theme(legend.position = "right")
CombinePlots(plots = list(plot1, plot2))

# Subset: remove potential empties, multiplets and broken cells based on parameters
seurat_obj <- PercentageFeatureSet(seurat_obj, "^mt-", col.name = "percent_mito")
seurat_obj= seurat_obj[, seurat_obj$nFeature_Spatial > mingenes & seurat_obj$percent_mito < mitocutoff]

#Sctransform normalizes the data, detects high-variance features, and stores the data in the SCT assay.
seurat_obj <- SCTransform(seurat_obj, assay = "Spatial", verbose = FALSE)

# close the pdf
dev.off() 

# Save the Robj for the next tool
save(seurat_obj, file="seurat_obj_sctransform.Robj")


## EOF
