# TOOL single-cell-seurat-filter-regress.R: "Seurat -Filtering, regression and detection of variable genes" (This tool filters out cells, normalises the data and regresses out uninteresting sources of variation in gene expression. It then detects highly variable genes across the single cells. You can use the plots from the Setup tool to estimate the parameter values.) 
# INPUT OPTIONAL seurat_obj.Robj: "Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL Dispersion_plot.pdf 
# OUTPUT OPTIONAL seurat_obj_2.Robj
# PARAMETER OPTIONAL mingenes: "Keep cells which express at least this many genes" TYPE INTEGER DEFAULT 200 (The cells need to have expressed at least this many genes.)
# PARAMETER OPTIONAL genecountcutoff: "Filter out cells which have higher unique gene count" TYPE INTEGER DEFAULT 2500 (Filter out potential multiplets, that is, cells that have more than this many unique gene counts.)
# PARAMETER OPTIONAL mitocutoff: "Filter out cells which have higher mitochondrial transcript percentage" TYPE DECIMAL FROM 0 TO 100 DEFAULT 5 (Filter out cells where the percentage of mitochondrial transcripts is higher than this.)
# PARAMETER OPTIONAL num.features: "Number of variable features to return" TYPE INTEGER DEFAULT 2000 (Number of features to select as top variable features, i.e. how many features returned.)
# PARAMETER OPTIONAL lognorm: "Perform log normalization" TYPE [T:yes, F:no] DEFAULT T (Select NO only if your data is already log transformed. For raw data, select YES.)
# PARAMETER OPTIONAL totalexpr: "Scale factor in the log normalization" TYPE INTEGER DEFAULT 10000 (Scale each cell to this total number of molecules before log normalization. Used in normalisation step.)
# PARAMETER OPTIONAL filter.cell.cycle: "Filter out cell cycle differences" TYPE [no:no, all.diff:"all differences", diff.phases:"the difference between the G2M and S phase scores"] DEFAULT no (Choose to remove all signal associated with cell cycle, or the difference between the G2M and S phase scores. More info in the manual page under Help. )
# RUNTIME R-3.4.3


# 2017-06-06 ML
# 2017-07-05 ML split into separate tool
# 2018-01-11 ML update Seurat version to 2.2.0
# 2018-02-19 ML Cell cycle filtering, plot dispersion plots without gene names and both scaled & non-scaled version, add parameters for normalisation
# 2018-11-07 ML Update the cutoffs for variable genes: xmin = 0.0125->0.1, x.high.cutoff = 3 -> 8, y.cutoff = 0.5 -> 1
# 2018-12-20 ML ScaleData: scale all changes in one command
# 2019-03-13 EK Removed names of variable genes from dispersion plot
# 2019-05-22 ML update Seurat version to 3.0


library(Seurat)
library(dplyr)
library(Matrix)
library(gplots)

# Load the R-Seurat-object (called seurat_obj)
load("seurat_obj.Robj")

# Open the pdf file for plotting
pdf(file="Dispersion_plot.pdf", width=13, height=7) 

# For cell cycle filtering, read in a list of cell cycle markers, from Tirosh et al, 2015
cc.genes <- readLines(con = file.path(chipster.tools.path, "seurat/regev_lab_cell_cycle_genes.txt"))
# We can segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]

# Subset
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > mingenes & nFeature_RNA < genecountcutoff & percent.mt < mitocutoff)	

if (lognorm=="T") {
	seurat_obj <- NormalizeData(object = seurat_obj, normalization.method = "LogNormalize", 
			scale.factor = totalexpr) 
}
# Detection of variable genes across the single cells
# FindVariableFeatures function identifies features that are outliers on a 'mean variability plot'.
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = num.features)

# Dispersion plot:
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seurat_obj), 10)

# Plot variable features with and without labels
plot1 <- VariableFeaturePlot(seurat_obj)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))

textplot(paste("\v \v Number of \n \v \v variable \n \v \v genes: \n \v \v", length(VariableFeatures(object = seurat_obj)), " \n  \n \v \v Number of \n \v \v cells: \n \v \v", length(colnames(x = seurat_obj))), halign="center", valign="center", cex=2)

## Scaling the data & cell cycle filtering:

# if filter.cell.cycle == "no", also anyhow RunPCA needs to access Object@scale.data:
seurat_obj <- ScaleData(object = seurat_obj, vars.to.regress = c("nCount_RNA", "percent.mt"), display.progress = FALSE)

if( filter.cell.cycle != "no" ) {	
	# Cell cycle genes, get the scores & visualize:
	# Note: in the very beginning we read in the table and set s.genes and gm2.genes
	# http://satijalab.org/seurat/cell_cycle_vignette.html#regress-out-cell-cycle-scores-during-data-scaling

	# Check that there were some S or G2M genes in the list of variable genes:
	if (length(s.genes[!is.na(match(s.genes, VariableFeatures(object = seurat_obj)))]) <1 && length(g2m.genes[!is.na(match(g2m.genes, VariableFeatures(object = seurat_obj)))]) <1 ) {
		stop(paste('CHIPSTER-NOTE: ', "There were no enough cell cycle genes for correction in the list of variable genes."))
	
	} else{
		seurat_obj <- CellCycleScoring(object = seurat_obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
		# Visualize the distribution of cell cycle markers across
		# Visualize in PCA:
		# PCA plot 1: before filtering cell cycle effect
		seurat_obj <- RunPCA(seurat_obj, features = c(s.genes, g2m.genes))
		plot1 <- DimPlot(seurat_obj, plot.title = "PCA on cell cycle genes") 
		# PCAPlot(seurat_obj, plot.title = "PCA on cell cycle genes") 
		# Remove the cell cycle scores:

		# Option 1: remove all the difference:
		if (filter.cell.cycle == "all.diff"){
			seurat_obj <- ScaleData(object = seurat_obj, vars.to.regress = c("S.Score", "G2M.Score", "nCount_RNA", "percent.mt"), 
					display.progress = FALSE)	
			# PCA plot 2A: after filtering, all:	
			seurat_obj <- RunPCA(object = seurat_obj, features = c(s.genes, g2m.genes), do.print = FALSE)
			plot2 <- DimPlot(seurat_obj, plot.title = "After cell cycle correction (method: remove all)")
			# PCAPlot(seurat_obj, plot.title = "After cell cycle correction (method: remove all)")

		# Option 2: regressing out the difference between the G2M and S phase scores:	
		}else if (filter.cell.cycle == "diff.phases"){
			seurat_obj$CC.Difference <- seurat_obj$S.Score - seurat_obj$G2M.Score
			seurat_obj <- ScaleData(object = seurat_obj, vars.to.regress = c("CC.Difference", "nCount_RNA", "percent.mt"), display.progress = FALSE)			
			# PCA plot 2B: after filtering, difference:	
			seurat_obj <- RunPCA(object = seurat_obj, features  = c(s.genes, g2m.genes), do.print = FALSE)
			plot2 <- DimPlot(seurat_obj, plot.title = "After cell cycle correction (method: difference between G2M and S phases)")
			# PCAPlot(seurat_obj, plot.title = "After cell cycle correction (method: difference between G2M and S phases)")
	
		}
		CombinePlots(plots = list(plot1, plot2))
	} 
} 
dev.off() # close the pdf

# Save the Robj for the next tool
save(seurat_obj, file="seurat_obj_2.Robj")

## EOF
