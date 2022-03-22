# TOOL single-cell-seurat-split-dot-plot.R: "Seurat v4 -Visualize genes with cell type specific responses in two samples" (This tool gives you plots showing user defined markers/genes across the conditions. This tool can be used for two sample combined Seurat objects.) 
# INPUT OPTIONAL combined_seurat_obj.Robj: "Combined Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL split_dot_plot.pdf
# PARAMETER markers: "Markers to plot" TYPE STRING DEFAULT "CD3D, CREM, HSPH1, SELL, GIMAP5" (Name of the marker genes you wish to plot, separated by comma. Please note that the gene names here are case sensitive, so check from your gene lists how the names are typed, e.g. CD3D vs Cd3d.)
# PARAMETER OPTIONAL reduction.method: "Visualisation with tSNE, UMAP or PCA" TYPE [umap:UMAP, tsne:tSNE, pca:PCA] DEFAULT umap (Which dimensionality reduction to use.)
# PARAMETER OPTIONAL plotting.order.used: "Plotting order of cells based on expression" TYPE [TRUE:yes, FALSE:no] DEFAULT FALSE (Plot cells in the the order of expression. Can be useful to turn this on if cells expressing given feature are getting buried.)
# IMAGE comp-20.04-r-deps
# RUNTIME R-4.1.0-single-cell


# 2018-16-05 ML
# 09.07.2019 ML Seurat v3
# 2020-06-17 ML Add ridge plot visualisation
# 2021-10-04 ML Update to Seurat v4

# For testing (not run):
# markers.to.plot <- c("CD3D", "CREM", "HSPH1", "SELL", "GIMAP5", "CACYBP", "GNLY", 
#		"NKG7", "CCL5", "CD8A", "MS4A1", "CD79A", "MIR155HG", "NME1", "FCGR3A", 
#		"VMO1", "CCL2", "S100A9", "HLA-DQA1", "GPR183", "PPBP", "GNG11", "HBA2", 
#		"HBB", "TSPAN13", "IL3RA", "IGJ")

# for UMAP:
library(reticulate)
use_python("/opt/chipster/tools/miniconda3/envs/chipster_tools/bin/python")

library(Seurat)

# Load the R-Seurat-object:
load("combined_seurat_obj.Robj")

if (exists("seurat_obj")) {
  data.combined <- seurat_obj
}


DefaultAssay(data.combined) <- "RNA" # this is very crucial.
# Store idents:
stored_idents <- Idents(data.combined)

markers.to.plot <- unlist(strsplit(markers, ", "))

# Sanity check: are the requested genes available in the data:
all.genes <- rownames(x = data.combined)
match(markers.to.plot, all.genes)
# if one of the genes is not in the list, print error message:
if (!all(!is.na(match(markers.to.plot, all.genes)))) { 
  not.found <- markers.to.plot[is.na(match(markers.to.plot, all.genes))==TRUE]
  stop(paste('CHIPSTER-NOTE: ', "The gene you requested was not found in this dataset:", not.found, " "))
  }

# pdf(file="split_dot_plot.pdf", , width=13, height=7)  # open pdf
pdf(file="split_dot_plot.pdf", width=12, height=12)  # open pdf


# Dot plot:
# DotPlot(data.combined, features = rev(markers.to.plot), cols = c("blue", "red"), dot.scale = 8, split.by = "stim") + RotatedAxis()
# Check how many samples there are and choose as many colors:
number.of.samples <- length(levels(as.factor((data.combined$stim))))
colors.for.samples <- rainbow(number.of.samples)
DotPlot(data.combined, features = rev(markers.to.plot), cols = colors.for.samples, dot.scale = 8, split.by = "stim") + RotatedAxis()


# Feature plot:
# Show in which cluster the genes are active
FeaturePlot(data.combined, features = markers.to.plot, min.cutoff = "q9", reduction=reduction.method, order=as.logical(plotting.order.used)) 

# Compare between the treatments:
# These plots get squeezed when there are many samples, and are at some point very difficult to read.
# That is why we subset the object so that only 4 samples are shown on each page.


number.of.samples <- length(levels(as.factor((data.combined$stim))))
sample.names <- levels(as.factor((data.combined$stim)))
Idents(data.combined) <- "stim"

if (number.of.samples <= 4){
FeaturePlot(data.combined, features = markers.to.plot, split.by = "stim", max.cutoff = 3, cols = c("grey", "blue"), reduction=reduction.method, order=as.logical(plotting.order.used))
  }
# If there are more than 4 samples, lets split them in multiple pages, using subsetting.
if (number.of.samples > 4) {
  # Using i as a bookmark, which sample we are now working with.
  i <- 1
  # Looping and printing pages until all samples are printed (rounding always up, last page can have only 1 sample.)
  for (j in 1:ceiling(number.of.samples/4) ) { 
    # If i is less than # of samples we have, keep subsetting and printing.
    if (i < number.of.samples) { 
      subset.of.samples <- subset(data.combined, idents = as.vector(sample.names[i:(i+3)] ) )
      Idents(subset.of.samples) <- "stim"
      # Need to save and print the plots for them to actually go to pdf:
      feat.plot <- FeaturePlot(subset.of.samples, features = markers.to.plot, split.by = "stim", max.cutoff = 3, cols = c("grey", "blue"), reduction=reduction.method, order=as.logical(plotting.order.used))
      print(feat.plot)
      i <- i+4
    }
  }
}

# FeatureHeatmap(data.combined, features.plot = markers.to.plot, group.by = "stim", pt.size = 0.25, key.position = "top", 
#		max.exp = 3)

DefaultAssay(data.combined) <- "RNA" # this is very crucial.
Idents(data.combined) <- stored_idents # Return the original idents

## Comparison violin plot:
data.combined$celltype <- Idents(data.combined)
plots <- VlnPlot(data.combined, features = markers.to.plot, split.by = "stim", group.by = "celltype",  pt.size = 0, combine = FALSE)
CombinePlots(plots = plots, ncol = 1)

## Ridge plot:
RidgePlot(data.combined, features = markers.to.plot, ncol = 2)

dev.off() # close the pdf

## EOF



