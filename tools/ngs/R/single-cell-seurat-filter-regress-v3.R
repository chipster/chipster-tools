# TOOL single-cell-seurat-filter-regress-v3.R: "Seurat v3 -Filter cells, normalize, regress and detect variable genes" (This tool filters out dead cells, empties and doublets. It then normalizes gene expression values and detects highly variable genes across the cells. Finally, it scales the data and regresses out unwanted variation based on the number of UMIs and mitochondrial transcript percentage. You can also choose to regress out variation due to cell cycle heterogeneity.) 
# INPUT OPTIONAL seurat_obj.Robj: "Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL Dispersion_plot.pdf 
# OUTPUT OPTIONAL seurat_obj_preprocess.Robj
# OUTPUT OPTIONAL log.txt
# PARAMETER OPTIONAL mingenes: "Filter out cells which have less than this many genes expressed" TYPE INTEGER DEFAULT 200 (Filter out empties. The cells to be kept must express at least this number of genes.)
# PARAMETER OPTIONAL genecountcutoff: "Filter out cells which have more than this many genes expressed" TYPE INTEGER DEFAULT 2500 (Filter out multiplets. The cells to be kept must express less than this number of genes.)
# PARAMETER OPTIONAL mitocutoff: "Filter out cells which have higher mitochondrial transcript percentage" TYPE DECIMAL FROM 0 TO 100 DEFAULT 5 (Filter out dead cells. The cells to be kept must have lower percentage of mitochondrial transcripts than this.)
# PARAMETER OPTIONAL lognorm: "Perform global scaling normalization" TYPE [T:yes, F:no] DEFAULT T (For raw data, select yes.)
# PARAMETER OPTIONAL totalexpr: "Scaling factor in the normalization" TYPE INTEGER DEFAULT 10000 (Scale each cell to this total number of transcripts.)
# PARAMETER OPTIONAL num.features: "Number of variable genes to return" TYPE INTEGER DEFAULT 2000 (Number of features to select as top variable features, i.e. how many features returned.)
# PARAMETER OPTIONAL filter.cell.cycle: "Regress out cell cycle differences" TYPE [no:no, all.diff:"all differences", diff.phases:"the difference between the G2M and S phase scores"] DEFAULT no (Would you like to regress out cell cycle scores during data scaling? If yes, should all signal associated with cell cycle be removed, or only the difference between the G2M and S phase scores.)
# IMAGE comp-20.04-r-base
# RUNTIME R-4.1.0-single-cell


# RUNTIME R-3.6.1-single-cell


# 2017-06-06 ML
# 2017-07-05 ML split into separate tool
# 2018-01-11 ML update Seurat version to 2.2.0
# 2018-02-19 ML Cell cycle filtering, plot dispersion plots without gene names and both scaled & non-scaled version, add parameters for normalisation
# 2018-11-07 ML Update the cutoffs for variable genes: xmin = 0.0125->0.1, x.high.cutoff = 3 -> 8, y.cutoff = 0.5 -> 1
# 2018-12-20 ML ScaleData: scale all changes in one command
# 2019-03-13 EK Removed names of variable genes from dispersion plot
# 2019-06-13 EK Add normalization in tool name, changed the order of parameters
# 2019-05-22 ML Update Seurat version to 3.0
# 2020-06-22 ML Update description
# 2020-07-02 ML Always compute the cell-cycle scoring and plot the PCA
# 2020-10-10 EK Update name, description and parameter order
# 2020-07-02 ML Remove the plot titles, as they started giving errors. Fix the IF structure.

library(Seurat)
library(dplyr)
library(Matrix)
library(gplots)
library(ggplot2)

# Load the R-Seurat-object (called seurat_obj)
load("seurat_obj.Robj")

# Open the pdf file for plotting
pdf(file="Dispersion_plot.pdf", width=13, height=7) 

# For cell cycle filtering, read in a list of cell cycle markers, from Tirosh et al, 2015
cc.genes <- readLines(con = file.path(chipster.tools.path, "seurat/regev_lab_cell_cycle_genes.txt"))
# We can segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]

# Subset: remove potential empties, multiplets and broken cells based on parameters.
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > mingenes & nFeature_RNA < genecountcutoff & percent.mt < mitocutoff)	

# Normalisation, scaling & finding variables genes:

if (lognorm=="T") {
	seurat_obj <- NormalizeData(object = seurat_obj, normalization.method = "LogNormalize", scale.factor = totalexpr) 
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

## Scaling: 
seurat_obj <- ScaleData(object = seurat_obj, vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = FALSE)

# Cell cycle stage scoring & PCA plot:
# Note: in the very beginning we read in the table and set s.genes and gm2.genes
# http://satijalab.org/seurat/cell_cycle_vignette.html#regress-out-cell-cycle-scores-during-data-scaling

# Check that there were some S or G2M genes in the list of variable genes:
if (length(s.genes[!is.na(match(s.genes, VariableFeatures(object = seurat_obj)))]) <1 && length(g2m.genes[!is.na(match(g2m.genes, VariableFeatures(object = seurat_obj)))]) <1 ) {
	# stop(paste('CHIPSTER-NOTE: ', "There were not enough cell cycle genes for correction in the list of variable genes."))
	# Write a log file (instead of ending with a Chipster-note, because we want the tool to finish and give the plot, when possible.)
	fileConn<-file("log.txt")
	writeLines(c("There are not enough cell cycle genes for correction in the list of variable genes."), fileConn)
	close(fileConn)
	
} else{
	seurat_obj <- CellCycleScoring(object = seurat_obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

	# Visualize in PCA:
	# PCA plot 1: without/before filtering cell cycle effect
	seurat_obj <- RunPCA(seurat_obj, features = c(s.genes, g2m.genes))
	plot1 <- DimPlot(seurat_obj) #, plot.title = "PCA on cell cycle genes")   reduction = pca

	# Cell cycle stage filtering:
	if( filter.cell.cycle != "no" ) {	
	# Remove the cell cycle scores:

		# Option 1: remove all the difference:
		if (filter.cell.cycle == "all.diff"){
			# Scale again, but this time including also the cell cycle scores
			seurat_obj <- ScaleData(object = seurat_obj, vars.to.regress = c("S.Score", "G2M.Score", "nCount_RNA", "percent.mt"), verbose = FALSE)

			# PCA plot 2A: after filtering, all:	
			seurat_obj <- RunPCA(object = seurat_obj, features = c(s.genes, g2m.genes), do.print = FALSE)
			ggtitle("After cell cycle correction (method: remove all)")
			plot2 <- DimPlot(seurat_obj)#, plot.title = "After cell cycle correction (method: remove all)")
			CombinePlots(plots = list(plot1, plot2))
		# Option 2: regressing out the difference between the G2M and S phase scores:	
		}else if (filter.cell.cycle == "diff.phases"){
			seurat_obj$CC.Difference <- seurat_obj$S.Score - seurat_obj$G2M.Score
			seurat_obj <- ScaleData(object = seurat_obj, vars.to.regress = c("CC.Difference", "nCount_RNA", "percent.mt"), verbose = FALSE)			
			# PCA plot 2B: after filtering, difference:	
			seurat_obj <- RunPCA(object = seurat_obj, features  = c(s.genes, g2m.genes), do.print = FALSE)
			ggtitle("After cell cycle correction (method: difference between G2M and S phases)")
			plot2 <- DimPlot(seurat_obj) #, plot.title = "After cell cycle correction (method: difference between G2M and S phases)")	
			CombinePlots(plots = list(plot1, plot2))
		}
	# just plot the 1 PCA plot, if no filtering:
	} else { 
		DimPlot(seurat_obj) #, plot.title = "PCA on cell cycle genes") 
	} 
} 

dev.off() # close the pdf

# Save the Robj for the next tool
save(seurat_obj, file="seurat_obj_preprocess.Robj")

## EOF
