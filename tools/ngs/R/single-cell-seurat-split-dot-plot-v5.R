# TOOL single-cell-seurat-split-dot-plot-v5.R: "Seurat v5 -Visualize genes with cell type specific responses in multiple samples" (This tool gives you plots showing user defined markers/genes across the conditions. This tool can be used for Seurat objects containing two or more samples.)
# INPUT OPTIONAL combined_seurat_obj.Robj: "Combined Seurat object" TYPE GENERIC
# INPUT OPTIONAL markers.txt: "Optional text file of the markers to plot" TYPE GENERIC (The names of the marker genes you wish to plot can also be given in the form of a text file, separated by comma. Please note that the gene names here are case sensitive, so check from your gene lists how the names are typed, e.g. CD3D vs Cd3d. In case the text file is provided, the markers to plot parameter is ignored.)
# OUTPUT OPTIONAL split_dot_plot.pdf
# PARAMETER OPTIONAL normalisation.method: "Normalisation method used previously" TYPE [LogNormalize:"Global scaling normalization", SCT:"SCTransform"] DEFAULT LogNormalize (Which normalisation method was used in preprocessing, Global scaling normalization \(default, NormalizeData function used\) or SCTransform.)
# PARAMETER OPTIONAL markers: "Markers to plot" TYPE STRING DEFAULT "CD3D, CREM, HSPH1, SELL, GIMAP5" (Name of the marker genes you wish to plot, separated by comma. Please note that the gene names here are case sensitive, so check from your gene lists how the names are typed, e.g. CD3D vs Cd3d.)
# PARAMETER OPTIONAL scale.data: "Scale data in split dot plot" TYPE [TRUE:yes, FALSE:no] DEFAULT FALSE (Determine whether the data is scaled in the split dot plot. By default, the raw expression data is used.)
# PARAMETER OPTIONAL reduction.method: "Visualisation with tSNE, UMAP or PCA in feature plot" TYPE [umap:UMAP, tsne:tSNE, pca:PCA] DEFAULT umap (Which dimensionality reduction to use in the feature plot.)
# PARAMETER OPTIONAL plotting.order.used: "Plotting order of cells based on expression in feature plot" TYPE [TRUE:yes, FALSE:no] DEFAULT FALSE (Plot cells in the the order of expression. Can be useful to turn this on if cells expressing a given feature are getting buried.)
# PARAMETER OPTIONAL color.scale: "Determine color scale based on all features in feature plot" TYPE [all:yes, feature:no] DEFAULT feature (Determine whether the color scales in the feature plots are based on all genes or individual genes. By default, the color scale is determined for each gene individually and may differ between genes. If you wish to compare gene expression between different genes, it is useful to set this parameter to "yes" so that the color scale is the same for all genes.)
# RUNTIME R-4.3.2-single-cell
# SLOTS 1


# 2018-16-05 ML
# 09.07.2019 ML Seurat v3
# 2020-06-17 ML Add ridge plot visualisation
# 2021-10-04 ML Update to Seurat v4
# 2023-12-12 IH Seurat v5

# For testing (not run):
# markers.to.plot <- c("CD3D", "CREM", "HSPH1", "SELL", "GIMAP5", "CACYBP", "GNLY",
# 		"NKG7", "CCL5", "CD8A", "MS4A1", "CD79A", "MIR155HG", "NME1", "FCGR3A",
# 		"VMO1", "CCL2", "S100A9", "HLA-DQA1", "GPR183", "PPBP", "GNG11", "HBA2",
# 		"HBB", "TSPAN13", "IL3RA", "IGJ")

# for UMAP:
# library(reticulate)
# Sys.setenv(RETICULATE_PYTHON = paste(chipster.tools.path, "/miniconda3/envs/chipster_tools/bin/python"))
# use_python("/opt/chipster/tools/miniconda3/envs/chipster_tools/bin/python")

library(Seurat)
library(readr)
options(Seurat.object.assay.version = "v5")

# for the fileOk function
source(file.path(chipster.common.lib.path, "tool-utils.R"))

# Load the R-Seurat-object:
load("combined_seurat_obj.Robj")

if (exists("seurat_obj")) {
  data.combined <- seurat_obj
}

# Use markers text file if provided, else the markers to plot parameter is used
if (fileOk("markers.txt", 0)) {
  markers <- read_file("markers.txt")
  markers.to.plot <- trimws(unlist(strsplit(markers, ",")))
} else {
  markers.to.plot <- trimws(unlist(strsplit(markers, ",")))
}

# In case some other type of assay is set:
if (normalisation.method == "SCT") {
  DefaultAssay(data.combined) <- "SCT"
} else {
  DefaultAssay(data.combined) <- "RNA"
}
# Store idents:
stored_idents <- Idents(data.combined)

# Sanity check: are the requested genes available in the data:
all.genes <- rownames(x = data.combined)
match(markers.to.plot, all.genes)
# if more than one of the genes is not in the list, print error message:
if (sum(is.na((match(markers.to.plot, all.genes)))) > 1) {
  not.found <- (markers.to.plot[is.na(match(markers.to.plot, all.genes)) == TRUE])
  not.found <- paste(not.found, collapse = ",")
  stop(paste("CHIPSTER-NOTE: ", "The genes you requested were not found in this dataset:", not.found))
}

# continue even if one gene is missing
if (!all(!is.na(match(markers.to.plot, all.genes)))) {
  not.found <- markers.to.plot[is.na(match(markers.to.plot, all.genes)) == TRUE]
  print(paste("Continuing the visualization without the one gene not found: ", not.found))
  markers.to.plot <- markers.to.plot[!is.na(match(markers.to.plot, all.genes))]
}

# pdf(file="split_dot_plot.pdf", , width=13, height=7)  # open pdf
pdf(file = "split_dot_plot.pdf", width = 12, height = 12) # open pdf


# Dot plot:
# DotPlot(data.combined, features = rev(markers.to.plot), cols = c("blue", "red"), dot.scale = 8, split.by = "stim") + RotatedAxis()
# Check how many samples there are and choose as many colors:
number.of.samples <- length(levels(as.factor((data.combined$stim))))
colors.for.samples <- rainbow(number.of.samples)
DotPlot(data.combined, features = rev(markers.to.plot), cols = colors.for.samples, dot.scale = 8, split.by = "stim", scale = as.logical(scale.data)) + RotatedAxis()


# Feature plot:
# Show in which cluster the genes are active
FeaturePlot(data.combined, features = markers.to.plot, min.cutoff = "q9", reduction = reduction.method, order = as.logical(plotting.order.used), keep.scale = color.scale)

# Compare between the treatments:
# These plots get squeezed when there are many samples, and are at some point very difficult to read.
# That is why we subset the object so that only 4 samples are shown on each page.


number.of.samples <- length(levels(as.factor((data.combined$stim))))
sample.names <- levels(as.factor((data.combined$stim)))
Idents(data.combined) <- "stim"

if (number.of.samples <= 4) {
  FeaturePlot(data.combined, features = markers.to.plot, split.by = "stim", max.cutoff = 3, cols = c("grey", "blue"), reduction = reduction.method, order = as.logical(plotting.order.used), keep.scale = color.scale)
}
# If there are more than 4 samples, lets split them in multiple pages, using subsetting.
if (number.of.samples > 4) {
  # Using i as a bookmark, which sample we are now working with.
  i <- 1
  # Looping and printing pages until all samples are printed (rounding always up, last page can have only 1 sample.)
  for (j in 1:ceiling(number.of.samples / 4)) {
    # If i is less than # of samples we have, keep subsetting and printing.
    if (i < number.of.samples) {
      samples.of.this.round <- as.vector(sample.names[i:(i + 3)])
      samples.of.this.round <- samples.of.this.round[!is.na(samples.of.this.round)]
      subset.of.samples <- subset(data.combined, idents = samples.of.this.round)
      Idents(subset.of.samples) <- "stim"
      # Need to save and print the plots for them to actually go to pdf:
      feat.plot <- FeaturePlot(subset.of.samples, features = markers.to.plot, split.by = "stim", max.cutoff = 3, cols = c("grey", "blue"), reduction = reduction.method, order = as.logical(plotting.order.used), keep.scale = color.scale)
      print(feat.plot)
      i <- i + 4
    }
  }
}

# FeatureHeatmap(data.combined, features.plot = markers.to.plot, group.by = "stim", pt.size = 0.25, key.position = "top",
# 		max.exp = 3)
if (normalisation.method == "SCT") {
  DefaultAssay(data.combined) <- "SCT"
} else {
  DefaultAssay(data.combined) <- "RNA"
}
Idents(data.combined) <- stored_idents # Return the original idents

## Comparison violin plot:
data.combined$celltype <- Idents(data.combined)
# plots <- VlnPlot(data.combined, features = markers.to.plot, split.by = "stim", group.by = "celltype", pt.size = 0, combine = FALSE)
plots <- VlnPlot(data.combined, features = markers.to.plot, split.by = "stim", group.by = "celltype", pt.size = 0, combine = FALSE, split.plot = TRUE)
CombinePlots(plots = plots, ncol = 1)

## Ridge plot:
RidgePlot(data.combined, features = markers.to.plot, ncol = 2)

dev.off() # close the pdf

## EOF
