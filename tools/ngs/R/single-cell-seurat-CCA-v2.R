# TOOL single-cell-seurat-CCA-v2.R: "Seurat v2 -Combine two samples and perform CCA" (This tool can be used to do the canonical correlation analysis CCA and to combine two Seurat objects for later joined analysis. The two objects need to be named when created in the Seurat setup tool.) 
# INPUT OPTIONAL seurat_obj1.Robj: "Seurat object 1" TYPE GENERIC
# INPUT OPTIONAL seurat_obj2.Robj: "Seurat object 2" TYPE GENERIC
# OUTPUT OPTIONAL CCAplot.pdf
# OUTPUT OPTIONAL seurat_obj_combined.Robj
# PARAMETER OPTIONAL CCstocompute: "How many CCs to compute" TYPE INTEGER DEFAULT 20 (Number of canonical vectors to calculate)
# PARAMETER OPTIONAL CCstovisualise: "How many CCs to visualise as heatmaps" TYPE INTEGER DEFAULT 9 (How many canonical components you want to visualise as heatmaps.)
# RUNTIME R-3.4.3
# SLOTS 3


# SLOTS = 3: when testing at VM this tool required 18.8G)

# 2018-08-05 ML
# 2018-10-03 ML Add sample identifiers to cell barcodes (fix problem with same cell barcodes in two samples)

library(Seurat)

# Load the R-Seurat-objects (called seurat_obj -that's why we need to rename them here)
load("seurat_obj1.Robj")
seurat_obj1 <- seurat_obj
load("seurat_obj2.Robj")
seurat_obj2 <- seurat_obj

# Gene selection for input to CCA
group1 <- FindVariableGenes(seurat_obj1, do.plot = F)
group2 <- FindVariableGenes(seurat_obj2, do.plot = F)

g.1 <- head(rownames(group1@hvg.info), 1000)
g.2 <- head(rownames(group2@hvg.info), 1000)

genes.use <- unique(c(g.1, g.2))
genes.use <- intersect(genes.use, rownames(group1@scale.data))
genes.use <- intersect(genes.use, rownames(group2@scale.data))

# Perform CCA
# data.combined <- RunCCA(seurat_obj1, seurat_obj2, genes.use = genes.use, num.cc = 30)
# In case there happen to be same cell barcodes, add identifiers to barcode names:
# https://github.com/satijalab/seurat/issues/135
data.combined<-RunCCA(seurat_obj1, seurat_obj2, genes.use = genes.use, num.cc = CCstocompute,
		add.cell.id1="one", add.cell.id2="two")


# Visualize results of CCA plot CC1 versus CC2 and look at a violin plot
pdf(file="CCAplot.pdf", , width=13, height=7)  # open pdf
p1 <- DimPlot(object = data.combined, reduction.use = "cca", group.by = "stim", 
		pt.size = 0.5, do.return = TRUE)
p2 <- VlnPlot(object = data.combined, features.plot = "CC1", group.by = "stim", 
		do.return = TRUE)
plot_grid(p1, p2)

# MetageneBicorPlot examines a measure of correlation strength for each CC 
# and find that this statistic generally saturates after a reasonable number of CCs.
p3 <- MetageneBicorPlot(data.combined, grouping.var = "stim", dims.eval = 1:CCstocompute, 
		display.progress = FALSE)

# As with PC selection, it is often also useful to examine heatmaps of the top genes 
# driving each CC. 
DimHeatmap(object = data.combined, reduction.type = "cca", cells.use = 500, 
		dim.use = 1:CCstovisualise, do.balanced = TRUE)

dev.off() # close the pdf


# PrintDim(object = immune.combined, reduction.type = "cca", dims.print = 1:2, 
#		genes.print = 10)

# Save the Robj for the next tool
save(data.combined, file="seurat_obj_combined.Robj")

## EOF



