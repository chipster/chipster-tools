# TOOL single-cell-seurat-filter-regress.R: "BETA Seurat -Filtering, regression and detection of variable genes" (This tool filters out cells and regresses out uninteresting sources of variation in gene expression. It then detects highly variable genes across the single cells. PLEASE NOTE that you might need to run the tool couple of times, as setting the max and min limits to average expression and dispersion with the bottom three parameters is an iterative process. Start with some values, see how it goes and run the tool again with different parameters.) 
# INPUT OPTIONAL seurat_obj.Robj: "Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL Dispersion_plot.pdf 
# OUTPUT OPTIONAL seurat_obj_2.Robj
# PARAMETER OPTIONAL genecountcutoff: "Filter out cells which have higher unique gene count" TYPE INTEGER DEFAULT 2500 (Filter out potential multiplets, that is, cells that have more than this many unique gene counts.)
# PARAMETER OPTIONAL mitocutoff: "Filter out cells which have higher mitochondrial transcript ratio" TYPE DECIMAL FROM 0 TO 1 DEFAULT 0.05 (Filter out cells where the ratio of mitochondrial transcripts is higher than this.)
# PARAMETER OPTIONAL xlowcutoff: "Minimum average expression level for a variable gene, x min" TYPE DECIMAL DEFAULT 0.0125 (For limiting the selection of variable genes.)
# PARAMETER OPTIONAL xhighcutoff: "Maximum average expression level for a variable gene, x max" TYPE DECIMAL DEFAULT 3 (For limiting the selection of variable genes.)
# PARAMETER OPTIONAL ylowcutoff: "Minimum dispersion for a variable gene, y min" TYPE DECIMAL DEFAULT 0.5 (For limiting the selection of variable genes.)
# RUNTIME R-3.3.2


# PARAMETER OPTIONAL yhighcutoff: "Top cutoff on y-axis for identifying variable genes" TYPE DECIMAL DEFAULT Inf (For limiting the selection of variable genes.)
# OUTPUT OPTIONAL log.txt

# 2017-06-06 ML
# 2017-07-05 ML split into separate tool

library(Seurat)
library(dplyr)
library(Matrix)
library(gplots)


# Load the R-Seurat-object (called seurat_obj)
load("seurat_obj.Robj")


# Filter out cells that have unique gene counts over threshold (default = 2,500)
# Other type of filtering could be included as well.
seurat_obj <- SubsetData(seurat_obj, subset.name = "nGene", accept.high = genecountcutoff)
seurat_obj <- SubsetData(seurat_obj, subset.name = "percent.mito", accept.high = mitocutoff)

# Regress unwanted sources of variation
# Other type of regressing could be included as well: 
# for example, by batch (if applicable) & cell alignment rate (provided by DropSeq tools for Drop-seq data).
# Here: number of detected molecules and mito gene expression.
seurat_obj <- RegressOut(seurat_obj, latent.vars = c("nUMI", "percent.mito"))

# Detection of variable genes across the single cells
# Identifies genes that are outliers on a 'mean variability plot'. 
# First, uses a function to calculate average expression (fxn.x) and dispersion (fxn.y) for each gene. 
# Next, divides genes into num.bin (deafult 20) bins based on their average expression, 
# and calculates z-scores for dispersion within each bin. 
# The purpose of this is to identify variable genes while controlling for the strong relationship 
# between variability and average expression.
pdf(file="Dispersion_plot.pdf") # width=13, height=7
seurat_obj <- MeanVarPlot(seurat_obj ,fxn.x = expMean, fxn.y = logVarDivMean, x.low.cutoff = xlowcutoff, 
		x.high.cutoff = xhighcutoff, y.cutoff = ylowcutoff, do.contour = F)
# y.high.cutoff = yhighcutoff,
length(seurat_obj@var.genes)
textplot(paste("Number of variable genes:", length(seurat_obj@var.genes)))
dev.off() # close the pdf


# Save the Robj for the next tool
save(seurat_obj, file="seurat_obj_2.Robj")

## EOF

