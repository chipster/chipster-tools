# TOOL single-cell-seurat-split-dot-plot-v3.R: "Seurat v3 BETA -Visualize genes with cell type specific responses in two samples" (This tool gives you plots showing user defined markers/genes across the conditions. This tool can be used for two sample combined Seurat objects.) 
# INPUT OPTIONAL combined_seurat_obj.Robj: "Combined Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL split_dot_plot.pdf
# PARAMETER markers: "Markers to plot" TYPE STRING DEFAULT "CD3D, CREM, HSPH1, SELL, GIMAP5" (Name of the marker genes you wish to plot, separated by comma. Please note that the gene names here are case sensitive, so check from your gene lists how the names are typed, e.g. CD3D vs Cd3d.)
# PARAMETER OPTIONAL reduction.method: "Visualisation with tSNE, UMAP or PCA" TYPE [umap:UMAP, tsne:tSNE, pca:PCA] DEFAULT umap (Which dimensionality reduction to use.)
# RUNTIME R-3.6.1


# 2018-16-05 ML
# 09.07.2019 ML Seurat v3

# For testing (not run):
# markers.to.plot <- c("CD3D", "CREM", "HSPH1", "SELL", "GIMAP5", "CACYBP", "GNLY", 
#		"NKG7", "CCL5", "CD8A", "MS4A1", "CD79A", "MIR155HG", "NME1", "FCGR3A", 
#		"VMO1", "CCL2", "S100A9", "HLA-DQA1", "GPR183", "PPBP", "GNG11", "HBA2", 
#		"HBB", "TSPAN13", "IL3RA", "IGJ")

# for UMAP:
library(reticulate)
use_python("/opt/chipster/tools/miniconda3/envs/chipster_tools/bin/python")

library(Seurat)

# Load the R-Seurat-objects (called seurat_obj -that's why we need to rename them here)
load("combined_seurat_obj.Robj")
#combined_seurat_obj <- data.combined

markers.to.plot <- unlist(strsplit(markers, ", "))
# pdf(file="split_dot_plot.pdf", , width=13, height=7)  # open pdf
pdf(file="split_dot_plot.pdf", width=12, height=12)  # open pdf


# Dot plot:
DotPlot(data.combined, features = rev(markers.to.plot), cols = c("blue", "red"), dot.scale = 8, 
	split.by = "stim") + RotatedAxis()
	

# Feature plot:
# Show in which cluster the genes are active
FeaturePlot(data.combined, features = markers.to.plot, min.cutoff = "q9", reduction=reduction.method) 
# Compare between the treatments:
FeaturePlot(data.combined, features = markers.to.plot, split.by = "stim", max.cutoff = 3, cols = c("grey", "red"), reduction=reduction.method)

# FeatureHeatmap(data.combined, features.plot = markers.to.plot, group.by = "stim", pt.size = 0.25, key.position = "top", 
#		max.exp = 3)


## Comparison violin plot:
data.combined$celltype <- Idents(data.combined)
plots <- VlnPlot(data.combined, features = markers.to.plot, split.by = "stim", group.by = "celltype",  pt.size = 0, combine = FALSE)
CombinePlots(plots = plots, ncol = 1)

dev.off() # close the pdf

## EOF



