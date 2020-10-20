# TOOL single-cell-seurat-sctransform.R: "Seurat v3 -SCTransform: Filter cells, normalize, regress and detect variable genes" (This tool filters out dead cells, empties and doublets. It then normalizes gene expression values using the SCTransform method, detects highly variable genes, scales the data and regresses out unwanted variation based on the number of UMIs and mitochondrial transcript percentage. You can also choose to regress out variation due to cell cycle heterogeneity.)  
# INPUT OPTIONAL seurat_obj.Robj: "Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL Dispersion_plot.pdf 
# OUTPUT OPTIONAL seurat_obj_2.Robj
# PARAMETER OPTIONAL mingenes: "Filter out cells which have less than this many genes expressed" TYPE INTEGER DEFAULT 200 (Filter out empties. The cells to be kept must express at least this number of genes.)
# PARAMETER OPTIONAL genecountcutoff: "Filter out cells which have more than this many genes expressed" TYPE INTEGER DEFAULT 2500 (Filter out multiplets. The cells to be kept must express less than this number of genes.)
# PARAMETER OPTIONAL mitocutoff: "Filter out cells which have higher mitochondrial transcript percentage" TYPE DECIMAL FROM 0 TO 100 DEFAULT 5 (Filter out dead cells. The cells to be kept must have lower percentage of mitochondrial transcripts than this.)
# PARAMETER OPTIONAL num.features: "Number of variable genes to return" TYPE INTEGER DEFAULT 3000 (Number of features to select as top variable features, i.e. how many features returned. For SCTransform, the recommended default is 3000.)
# PARAMETER OPTIONAL filter.cell.cycle: "Regress out cell cycle differences" TYPE [no:no, all.diff:"all differences", diff.phases:"the difference between the G2M and S phase scores"] DEFAULT no (Would you like to regress out cell cycle scores during data scaling? If yes, should all signal associated with cell cycle be removed, or only the difference between the G2M and S phase scores.)
# RUNTIME R-3.6.1-single-cell


# 2020-06-17 ML
# 2020-10-11 EK Unified parameter descriptions with the corresponding normalization tool 

# Source: https://github.com/satijalab/seurat/issues/1679


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

# SCTransform:
# Note that this single command replaces NormalizeData, ScaleData, and FindVariableFeatures.

# Normalize data with SCTransform: 
seurat_obj <- SCTransform(seurat_obj, assay = 'RNA', new.assay.name = 'SCT', vars.to.regress = c('percent.mt', 'nFeature_RNA', 'nCount_RNA'), variable.features.n = num.features, verbose = FALSE)


# Dispersion plot:
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seurat_obj), 10)

# Plot variable features with and without labels
plot1 <- VariableFeaturePlot(seurat_obj)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))

textplot(paste("\v \v Number of \n \v \v variable \n \v \v genes: \n \v \v", length(VariableFeatures(object = seurat_obj)), " \n  \n \v \v Number of \n \v \v cells: \n \v \v", length(colnames(x = seurat_obj))), halign="center", valign="center", cex=2)

## Cell cycle filtering:

## if filter.cell.cycle == "no", also anyhow RunPCA needs to access Object@scale.data:
## seurat_obj <- ScaleData(object = seurat_obj, vars.to.regress = c("nCount_RNA", "percent.mt"), display.progress = FALSE)


if( filter.cell.cycle != "no" ) {	
	# Cell cycle genes, get the scores & visualize:
	# Note: in the very beginning we read in the table and set s.genes and gm2.genes
	# http://satijalab.org/seurat/cell_cycle_vignette.html#regress-out-cell-cycle-scores-during-data-scaling

	# Check that there were some S or G2M genes in the list of variable genes:
	if (length(s.genes[!is.na(match(s.genes, VariableFeatures(object = seurat_obj)))]) <1 && length(g2m.genes[!is.na(match(g2m.genes, VariableFeatures(object = seurat_obj)))]) <1 ) {
		stop(paste('CHIPSTER-NOTE: ', "There were no enough cell cycle genes for correction in the list of variable genes."))
	
	} else{

        # Perform cell cycle analysis (make sure to specify the "assay" parameter
		seurat_obj <- CellCycleScoring(object = seurat_obj, s.features = s.genes, g2m.features = g2m.genes, assay = 'SCT', set.ident = TRUE)
		
        # Visualize the distribution of cell cycle markers across
		# Visualize in PCA:
		# PCA plot 1: before filtering cell cycle effect
		seurat_obj <- RunPCA(seurat_obj, features = c(s.genes, g2m.genes))
		plot1 <- DimPlot(seurat_obj, plot.title = "PCA on cell cycle genes") 
		# PCAPlot(seurat_obj, plot.title = "PCA on cell cycle genes") 
		# Remove the cell cycle scores:

		# Option 1: remove all the difference:
		if (filter.cell.cycle == "all.diff"){
            # Normalise again with SCTransform, but this time including also the cell cycle scores
            seurat_obj <- SCTransform(seurat_obj, assay = 'RNA', new.assay.name = 'SCT', vars.to.regress = c('percent.mt', 'nFeature_RNA', 'nCount_RNA', 'S.Score', 'G2M.Score'), variable.features.n = num.features, verbose = FALSE)

            # PCA plot 2A: after filtering, all:	
			seurat_obj <- RunPCA(object = seurat_obj, features = c(s.genes, g2m.genes), do.print = FALSE)
			plot2 <- DimPlot(seurat_obj, plot.title = "After cell cycle correction (method: remove all)")
			# PCAPlot(seurat_obj, plot.title = "After cell cycle correction (method: remove all)")

		# Option 2: regressing out the difference between the G2M and S phase scores:	
		}else if (filter.cell.cycle == "diff.phases"){
			seurat_obj$CC.Difference <- seurat_obj$S.Score - seurat_obj$G2M.Score
            seurat_obj <- SCTransform(seurat_obj, assay = 'RNA', new.assay.name = 'SCT', vars.to.regress = c('percent.mt', 'nFeature_RNA', 'nCount_RNA', 'CC.Difference'), variable.features.n = num.features, verbose = FALSE)

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
