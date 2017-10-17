# TOOL single-cell-seurat-setup.R: "BETA Seurat -Setup and QC" (Setup the Seurat object and examine quality control plots. As an input, give either a tar package of a folder which contains the 10X output files, or a DGE matrix for DropSeq data. Please check that your input is assigned correctly under the parameters!)
# INPUT OPTIONAL files.tar: "tar package of 10X output files" TYPE GENERIC
# INPUT OPTIONAL dropseq.tsv: "DGE table from DropSeq" TYPE GENERIC
# OUTPUT OPTIONAL QCplots.pdf 
# OUTPUT OPTIONAL PCAplots.pdf
# OUTPUT OPTIONAL PCAgenes.txt
# OUTPUT OPTIONAL seurat_obj.Robj
# PARAMETER OPTIONAL project.name: "Project name for plotting" TYPE STRING DEFAULT Project_name (You can give your project a name. The name will appear on the plots.)
# PARAMETER OPTIONAL mincells: "Keep genes which are expressed in at least this many cells" TYPE INTEGER DEFAULT 3 (The genes need to be expressed in at least this many cells.)
# PARAMETER OPTIONAL mingenes: "Keep cells which express at least this many genes" TYPE INTEGER DEFAULT 200 (The cells need to have expressed at least this many genes.)
# PARAMETER OPTIONAL lognorm: "Perform log normalization" TYPE [T:yes, F:no] DEFAULT T (Select NO only if your data is already log transformed. For raw data, select YES.)
# PARAMETER OPTIONAL totalexpr: "Scale factor in the log normalization" TYPE INTEGER DEFAULT 10000 (Scale each cell to this total number of molecules before log normalization.)
# RUNTIME R-3.3.2


# PARAMETER OPTIONAL yhighcutoff: "Top cutoff on y-axis for identifying variable genes" TYPE DECIMAL DEFAULT Inf (For limiting the selection of variable genes.)
# OUTPUT OPTIONAL log.txt
# PARAMETER OPTIONAL genecountcutoff: "Unique gene counts per cell upper limit cutoff" TYPE INTEGER DEFAULT 2500 (Filter out potential multiplets, that is, cells that have more than this many unique gene counts.)
# PARAMETER OPTIONAL mitocutoff: "Mitochondrial genes percentage upper limit cutoff" TYPE DECIMAL FROM 0 TO 1 DEFAULT 0.05 (Filter out cells with higher than this percent of mitochondrial genes present.)


# 2017-06-06 ML
# 2017-07-05 ML split into separate tool

library(Seurat)
library(dplyr)
library(Matrix)
library(gplots)

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
seurat_obj <- new("seurat", raw.data = dat)
seurat_obj <- Setup(seurat_obj, min.cells = mincells, min.genes = mingenes, 
		do.logNormalize = lognorm, total.expr = totalexpr, project = project.name)


# QC
# % of mito genes
mito.genes <- grep("^MT-", rownames(seurat_obj@data), value = T)
percent.mito <- colSums(expm1(seurat_obj@data[mito.genes, ]))/colSums(expm1(seurat_obj@data))
seurat_obj <- AddMetaData(seurat_obj, percent.mito, "percent.mito")
# pdf plots
pdf(file="QCplots.pdf", , width=13, height=7) 
VlnPlot(seurat_obj, c("nGene", "nUMI", "percent.mito"), nCol = 3) 
par(mfrow = c(1, 2))
GenePlot(seurat_obj, "nUMI", "percent.mito")
# abline(h=mitocutoff, col="blue")
GenePlot(seurat_obj, "nUMI", "nGene")
# abline(h=genecountcutoff, col="blue")
#dev.off() # close the pdf


# Save the Robj for the next tool
save(seurat_obj, file="seurat_obj.Robj")

## EOF

