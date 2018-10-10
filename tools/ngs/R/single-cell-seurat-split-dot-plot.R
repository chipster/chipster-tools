# TOOL single-cell-seurat-split-dot-plot.R: "Seurat -Visualize genes with cell type specific responses in two samples" (This tool gives you plots showing user defined markers/genes across the conditions. This tool can be used for two sample combined Seurat objects.) 
# INPUT OPTIONAL combined_seurat_obj.Robj: "Combined Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL split_dot_plot.pdf
# PARAMETER markers: "Markers to plot" TYPE STRING DEFAULT "CD3D, CREM, HSPH1, SELL, GIMAP5" (Name of the marker genes you wish to plot, separated by comma.)
# RUNTIME R-3.4.3



# 2018-16-05 ML

# For testing (not run):
# markers.to.plot <- c("CD3D", "CREM", "HSPH1", "SELL", "GIMAP5", "CACYBP", "GNLY", 
#		"NKG7", "CCL5", "CD8A", "MS4A1", "CD79A", "MIR155HG", "NME1", "FCGR3A", 
#		"VMO1", "CCL2", "S100A9", "HLA-DQA1", "GPR183", "PPBP", "GNG11", "HBA2", 
#		"HBB", "TSPAN13", "IL3RA", "IGJ")

library(Seurat)

# Load the R-Seurat-objects (called seurat_obj -that's why we need to rename them here)
load("combined_seurat_obj.Robj")
#combined_seurat_obj <- data.combined

markers.to.plot <- unlist(strsplit(markers, ", "))
pdf(file="split_dot_plot.pdf", , width=13, height=7)  # open pdf
sdp <- SplitDotPlotGG(data.combined, grouping.var = "stim", genes.plot = rev(markers.to.plot), cols.use = c("blue", 
				"red"), x.lab.rot = T, plot.legend = T, dot.scale = 8, do.return = T)

FeatureHeatmap(data.combined, features.plot = markers.to.plot, group.by = "stim", pt.size = 0.25, key.position = "top", 
		max.exp = 3)

dev.off() # close the pdf

## EOF



