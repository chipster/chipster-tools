# TOOL single-cell-seurat-filter-regress.R: "BETA Seurat -Filtering, regression and detection of variable genes" (This tool filters out cells and regresses out uninteresting sources of variation in gene expression. It then detects highly variable genes across the single cells. PLEASE NOTE that you might need to run the tool couple of times, as setting the max and min limits to average expression and dispersion with the bottom three parameters is an iterative process. Start with some values, see how it goes and run the tool again with different parameters.) 
# INPUT OPTIONAL seurat_obj.Robj: "Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL Dispersion_plot.pdf 
# OUTPUT OPTIONAL seurat_obj_2.Robj
# PARAMETER OPTIONAL mingenes: "Keep cells which express at least this many genes" TYPE INTEGER DEFAULT 200 (The cells need to have expressed at least this many genes.)
# PARAMETER OPTIONAL genecountcutoff: "Filter out cells which have higher unique gene count" TYPE INTEGER DEFAULT 2500 (Filter out potential multiplets, that is, cells that have more than this many unique gene counts.)
# PARAMETER OPTIONAL mitocutoff: "Filter out cells which have higher mitochondrial transcript ratio" TYPE DECIMAL FROM 0 TO 1 DEFAULT 0.05 (Filter out cells where the ratio of mitochondrial transcripts is higher than this.)
# PARAMETER OPTIONAL xlowcutoff: "Minimum average expression level for a variable gene, x min" TYPE DECIMAL DEFAULT 0.0125 (For limiting the selection of variable genes.)
# PARAMETER OPTIONAL xhighcutoff: "Maximum average expression level for a variable gene, x max" TYPE DECIMAL DEFAULT 3 (For limiting the selection of variable genes.)
# PARAMETER OPTIONAL ylowcutoff: "Minimum dispersion for a variable gene, y min" TYPE DECIMAL DEFAULT 0.5 (For limiting the selection of variable genes.)
# PARAMETER OPTIONAL lognorm: "Perform log normalization" TYPE [T:yes, F:no] DEFAULT T (Select NO only if your data is already log transformed. For raw data, select YES.)
# PARAMETER OPTIONAL totalexpr: "Scale factor in the log normalization" TYPE INTEGER DEFAULT 10000 (Scale each cell to this total number of molecules before log normalization. Used in normalisation step.)
# PARAMETER OPTIONAL filter.cell.cycle: "Filter out cell cycle differences" TYPE [no:no, all.diff:"all differences", diff.phases:"the difference between the G2M and S phase scores"] DEFAULT no (Choose to remove all signal associated with cell cycle, or the difference between the G2M and S phase scores. More info in the manual page under Help. )
# RUNTIME R-3.4.3


# PARAMETER OPTIONAL yhighcutoff: "Top cutoff on y-axis for identifying variable genes" TYPE DECIMAL DEFAULT Inf (For limiting the selection of variable genes.)
# OUTPUT OPTIONAL log.txt

# 2017-06-06 ML
# 2017-07-05 ML split into separate tool
# 2018-01-11 ML update Seurat version to 2.2.0
# 2018-02-19 ML Cell cycle filtering, plot dispersion plots without gene names and both scaled & non-scaled version, add parameters for normalisation


library(Seurat)
library(dplyr)
library(Matrix)
library(gplots)


# Load the R-Seurat-object (called seurat_obj)
load("seurat_obj.Robj")

pdf(file="Dispersion_plot.pdf") # , width=13, height=7) 

# For cell cycle filtering, read in a list of cell cycle markers, from Tirosh et al, 2015
cc.genes <- readLines(con = file.path(chipster.tools.path, "seurat/regev_lab_cell_cycle_genes.txt"))
# We can segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]


# before or after cell cycle fixing?
seurat_obj <- FilterCells(object = seurat_obj, subset.names = c("nGene", "percent.mito"), 
		low.thresholds = c(mingenes, -Inf), high.thresholds = c(genecountcutoff, mitocutoff)) 

if (lognorm=="T") {
	seurat_obj <- NormalizeData(object = seurat_obj, normalization.method = "LogNormalize", 
			scale.factor = totalexpr) 
}
# Detection of variable genes across the single cells
# Identifies genes that are outliers on a 'mean variability plot'. 
# First, uses a function to calculate average expression (fxn.x) and dispersion (fxn.y) for each gene. 
# Next, divides genes into num.bin (deafult 20) bins based on their average expression, 
# and calculates z-scores for dispersion within each bin. 
# The purpose of this is to identify variable genes while controlling for the strong relationship 
# between variability and average expression.
seurat_obj <- FindVariableGenes(object = seurat_obj, mean.function = ExpMean, dispersion.function = LogVMR, 
		x.low.cutoff = xlowcutoff, x.high.cutoff = xhighcutoff, y.cutoff = ylowcutoff, do.text = FALSE, plot.both = TRUE )
# do.text	= Add text names of variable genes to plot (default is TRUE)
# plot.both	= Plot both the scaled and non-scaled graphs.
length(x = seurat_obj@var.genes)
seurat_obj <- ScaleData(object = seurat_obj, vars.to.regress = c("nUMI", "percent.mito"), display.progress = FALSE)
textplot(paste("\v \v Number of \n \v \v variable \n \v \v genes: \n \v \v", length(seurat_obj@var.genes)), halign="center", valign="center", cex=0.8)

if( filter.cell.cycle != "no" ) {
	
	# Cell cycle genes, get the scores & visualise:
	# Note: in the very beginning we read in the table and set s.genes and gm2.genes
	# http://satijalab.org/seurat/cell_cycle_vignette.html#regress-out-cell-cycle-scores-during-data-scaling

	# Check that there were some S or G2M genes in the list of variable genes:
	if (length(s.genes[!is.na(match(s.genes, seurat_obj@var.genes))]) <1 && length(g2m.genes[!is.na(match(g2m.genes, seurat_obj@var.genes))]) <1 ) {
		stop(paste('CHIPSTER-NOTE: ', "There were no enough cell cycle genes for correction in the list of variable genes."))
	} else{
		seurat_obj <- CellCycleScoring(object = seurat_obj, s.genes = s.genes, g2m.genes = g2m.genes, 
				set.ident = TRUE)
	#	 Visualize the distribution of cell cycle markers across
		# JoyPlot(object = seurat_obj, features.plot = c("PCNA", "TOP2A", "MCM6", "MKI67"), nCol = 2)  # Not all genes found, JoyPlot has been replaced with RidgePlot 
		# Visualize in PCA:
		seurat_obj <- RunPCA(object = seurat_obj, pc.genes = c(s.genes, g2m.genes), do.print = FALSE)
		PCAPlot(object = seurat_obj, plot.title = "PCA on cell cycle genes")

		# Remove the cell cycle scores:

		# Option 1: remove all the difference:
		if (filter.cell.cycle == "all.diff"){
			seurat_obj <- ScaleData(object = seurat_obj, vars.to.regress = c("S.Score", "G2M.Score"), 
					display.progress = FALSE)
			seurat_obj <- RunPCA(object = seurat_obj, pc.genes = c(s.genes, g2m.genes), do.print = FALSE)
			PCAPlot(object = seurat_obj, plot.title = "After cell cycle correction (method: remove all)") # HUOM, size
		# Option 2: regressing out the difference between the G2M and S phase scores:	
		}else if (filter.cell.cycle == "diff.phases"){
			seurat_obj@meta.data$CC.Difference <- seurat_obj@meta.data$S.Score - seurat_obj@meta.data$G2M.Score
			seurat_obj <- ScaleData(object = seurat_obj, vars.to.regress = "CC.Difference", display.progress = FALSE)
			seurat_obj <- RunPCA(object = seurat_obj, pc.genes = c(s.genes, g2m.genes), do.print = FALSE)
			PCAPlot(object = seurat_obj, plot.title = "After cell cycle correction (method: difference between G2M and S phases)") # HUOM size
		}
	} 
} 
dev.off() # close the pdf

# Save the Robj for the next tool
save(seurat_obj, file="seurat_obj_2.Robj")

## EOF
