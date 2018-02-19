# TOOL single-cell-seurat-setup.R: "BETA Seurat -Setup and QC" (Setup the Seurat object, quality control, filter and regress the cells, determine statistically significant principal components. As an input, give a .tar package of a folder which contains the 10X output files OR a DGE matrix for DropSeq data. Please check that your input is assigned correctly under the parameters!)
# INPUT OPTIONAL files.tar: "tar package of 10X output files" TYPE GENERIC
# INPUT OPTIONAL dropseq.tsv: "DGE table from DropSeq" TYPE GENERIC
# OUTPUT OPTIONAL QCplots.pdf 
# OUTPUT OPTIONAL PCAplots.pdf
# OUTPUT OPTIONAL PCAgenes.txt
# OUTPUT OPTIONAL seurat_obj.Robj
# PARAMETER OPTIONAL project.name: "Project name for plotting" TYPE STRING DEFAULT Project_name (You can give your project a name. The name will appear on the plots.)
# PARAMETER OPTIONAL mincells: "Include genes with detected expression in at least this many cells" TYPE INTEGER DEFAULT 3 (The genes need to be expressed in at least this many cells.)
# PARAMETER OPTIONAL mingenes: "Include cells where at least this many genes are detected" TYPE INTEGER DEFAULT 200 (The cells need to have expressed at least this many genes.)
# PARAMETER OPTIONAL genecountcutoff: "Unique gene counts per cell upper limit cutoff" TYPE INTEGER DEFAULT 2500 (Filter out potential mutliplets, that is, cells that have more than this many unique gene counts.)
# PARAMETER OPTIONAL mitocutoff: "Mitochondrial genes percentage upper limit cutoff" TYPE DECIMAL FROM 0 TO 1 DEFAULT 0.05 (Filter out cells with higher than this percent of mitochondrial genes present.)
# PARAMETER OPTIONAL xlowcutoff: "Bottom cutoff on x-axis for identifying variable genes" TYPE DECIMAL DEFAULT 0.0125 (For limiting the selection of variable genes.)
# PARAMETER OPTIONAL xhighcutoff: "Top cutoff on x-axis for identifying variable genes" TYPE DECIMAL DEFAULT 3 (For limiting the selection of variable genes.)
# PARAMETER OPTIONAL ylowcutoff: "Bottom cutoff on y-axis for identifying variable genes" TYPE DECIMAL DEFAULT 0.5 (For limiting the selection of variable genes.)
# PARAMETER OPTIONAL filter.cell.cycle: "Filter out cell cycle differences" TYPE [no:no, all.diff:"all differences", diff.phases:"the difference between the G2M and S phase scores"] DEFAULT no (Choose to remove all signal associated with cell cycle, or the difference between the G2M and S phase scores. More info in the manual page under Help. )
# RUNTIME R-3.4.3




## PARAMETER OPTIONAL yhighcutoff: "Top cutoff on y-axis for identifying variable genes" TYPE DECIMAL DEFAULT Inf (For limiting the selection of variable genes.)
## OUTPUT OPTIONAL log.txt

# MUISTA POISTAA KAKS TURHAKS JÄÄNYTTÄ PARAMETRIA!  do.logNormalize = lognorm, total.expr = totalexpr,
# PARAMETER OPTIONAL lognorm: "Perform log normalization" TYPE [T:yes, F:no] DEFAULT T (Select NO only if your data is already log transformed. For raw data, select YES.)
# PARAMETER OPTIONAL totalexpr: "Scale factor in the log normalization" TYPE INTEGER DEFAULT 10000 (Scale each cell to this total number of molecules before log normalization.)


# 2017-06-06 ML
# 2018-01-11 ML update Seurat version to 2.2.0
# 2018-01-22 ML add option to regress out cell cycle scores


library(Seurat)
library(dplyr)
library(Matrix)
library(gplots)

# For cell cycle filtering, read in a list of cell cycle markers, from Tirosh et al, 2015
cc.genes <- readLines(con = file.path(chipster.tools.path, "seurat/regev_lab_cell_cycle_genes.txt"))
# We can segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]

# If using DropSeq data:
if (file.exists("dropseq.tsv")){
	dat <- read.table("dropseq.tsv", header=T, sep="\t", row.names=1)
	
# If using 10X data:
}else if (file.exists("files.tar")){
	
	# Read the contents of the tar file into a list
	system("tar tf files.tar > tar.contents 2>> log.txt")
	file.list <- scan("tar.contents", what="", sep="\n")
	
	# Check that the input is a valid tar file
	if (length(file.list) == 0){
		stop(paste('CHIPSTER-NOTE: ', "It seems your input file is not a valid Tar package. Please check your input file."))
	}
	## Check if tar package contains folders
	#if (grepl("/", file.list[1])){
	#	stop(paste('CHIPSTER-NOTE: ', "It seems your Tar package contains folders. The input files need to be in the root of the package, not in subfolders."))
	#}
	
	# Open tar package
	system("tar xf files.tar 2>> log.txt")
	# Check the name of the directory 
	dir <- as.character(read.table("tar.contents")[1,1])
	
	# Load the data
	dat <- Read10X(dir)
}else{
	stop(paste('CHIPSTER-NOTE: ', "You need to provide either a 10X directory as a Tar package OR a DropSeq DGE as a tsv table. Please check your input file."))
}

# Initialize the Seurat object
# Huom, nämä siirtyi myöhemmäksi: do.logNormalize = lognorm, total.expr = totalexpr,

seurat_obj <- CreateSeuratObject(raw.data = dat, min.cells = mincells, min.genes = mingenes, 
		project = project.name)

# QC
# % of mito genes (note: they are named either "MT-CO1" or "mt-Co1", have to check both)
mito.genes <- grep(pattern ="^MT-", x = rownames(x =seurat_obj@data), value = T, ignore.case=T)
#percent.mito <- colSums(expm1(seurat_obj@data[mito.genes, ]))/colSums(expm1(seurat_obj@data))
percent.mito <- Matrix::colSums(seurat_obj@raw.data[mito.genes, ])/Matrix::colSums(seurat_obj@raw.data)
seurat_obj <- AddMetaData(object = seurat_obj, metadata = percent.mito, col.name = "percent.mito")
# pdf plots
pdf(file="QCplots.pdf", , width=13, height=7) 
VlnPlot(seurat_obj, c("nGene", "nUMI", "percent.mito"), nCol = 3) 
par(mfrow = c(1, 2))
GenePlot(seurat_obj, "nUMI", "percent.mito")
GenePlot(seurat_obj, "nUMI", "nGene")
#dev.off() # close the pdf



# tää ennen vai jälkeen cell cycle korjausta???
# HUOM, mingenes nyt kahdessa kohtaa...? pitäiskö olla erillinen parametri vai häh....? 200
seurat_obj <- FilterCells(object = seurat_obj, subset.names = c("nGene", "percent.mito"), 
		low.thresholds = c(mingenes, -Inf), high.thresholds = c(genecountcutoff, mitocutoff)) #Huom, parametrit muuttuu...???
seurat_obj <- NormalizeData(object = seurat_obj, normalization.method = "LogNormalize", 
		scale.factor = 10000) # olisko näille hyvä kuitenkin olla ne parametrit?
seurat_obj <- FindVariableGenes(object = seurat_obj, mean.function = ExpMean, dispersion.function = LogVMR, 
		x.low.cutoff = xlowcutoff, x.high.cutoff = xhighcutoff, y.cutoff = ylowcutoff)
length(x = seurat_obj@var.genes)
seurat_obj <- ScaleData(object = seurat_obj, vars.to.regress = c("nUMI", "percent.mito"), display.progress = FALSE)
textplot(paste("Number of variable \n genes: \n", length(seurat_obj@var.genes)), cex=0.8)


# Cell cycle genes, get the scores & visualise:
# Note: in the very beginning we read in the table and set s.genes and gm2.genes
# http://satijalab.org/seurat/cell_cycle_vignette.html#regress-out-cell-cycle-scores-during-data-scaling
seurat_obj <- CellCycleScoring(object = seurat_obj, s.genes = s.genes, g2m.genes = g2m.genes, 
		set.ident = TRUE)
# Visualize the distribution of cell cycle markers across
# JoyPlot(object = seurat_obj, features.plot = c("PCNA", "TOP2A", "MCM6", "MKI67"), nCol = 2)  # Not all genes found, JoyPlot has been replaced with RidgePlot 
# Visualize in PCA:
seurat_obj <- RunPCA(object = seurat_obj, pc.genes = c(s.genes, g2m.genes), do.print = FALSE)
PCAPlot(object = seurat_obj, plot.title = "PCA on cell cycle genes")

# Remove the cell cycle scores:
# Huom, lisää kysely, että löytyikö tarpeeksi geenejä listoilta! RunPCA herjaa:
#seurat_obj <- RunPCA(object = seurat_obj, pc.genes = c(s.genes, g2m.genes), do.print = FALSE)
#Error in irlba(A = t(x = data.use), nv = pcs.compute, ...) : 
#		max(nu, nv) must be positive
#Calls: RunPCA -> irlba
#Execution halted

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

dev.off() # close the pdf



# PCA
# The variable genes = genes in seurat_object@var.genes are used as input
seurat_obj <- RunPCA(object = seurat_obj, pc.genes = seurat_obj@var.genes, do.print = TRUE, pcs.print = 5, genes.print = 5)

#PCA genes in txt file
sink("PCAgenes.txt")
PrintPCA(seurat_obj, pcs.print = 1:5, genes.print = 5, use.full = TRUE)
sink()
# PDF
pdf(file="PCAplots.pdf", , width=9, height=12) 
VizPCA(seurat_obj, 1:2)
PCAPlot(seurat_obj, 1, 2)
PCHeatmap(seurat_obj, pc.use = 1, cells.use = 100, do.balanced = TRUE)
# fig.height=12,fig.width=9 
PCHeatmap(seurat_obj, pc.use = 1:12, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)
PCElbowPlot(seurat_obj)
dev.off() # close the pdf

# Save the Robj for the next tool
save(seurat_obj, file="seurat_obj.Robj")

## EOF

