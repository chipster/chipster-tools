# TOOL single-cell-seurat-setup-v3.R: "Seurat v3 -Setup and QC" (Setup the Seurat object, make quality control plots and filter out genes. There are two options for input files, please check that your input file is correctly assigned under the parameters. If you have 10X data, make a tar package containing the files genes.tsv, barcodes.tsv and matrix.mtx \(you can use the tool \"Utilities - Make a tar package\" for this\). Alternatively you can give a DGE matrix as input. If you are planning to combine samples later on, make sure you name them in this tool.)
# INPUT OPTIONAL files.tar: "tar package of 10X output files" TYPE GENERIC
# INPUT OPTIONAL dropseq.tsv: "DGE table" TYPE GENERIC
# OUTPUT OPTIONAL QCplots.pdf 
# OUTPUT OPTIONAL PCAplots.pdf
# OUTPUT OPTIONAL PCAgenes.txt
# OUTPUT OPTIONAL seurat_obj.Robj
# PARAMETER OPTIONAL project.name: "Project name for plotting" TYPE STRING DEFAULT Project_name (You can give your project a name. The name will appear on the plots. Do not use underscore _ in the names!)
# PARAMETER OPTIONAL mincells: "Keep genes which are expressed in at least this many cells" TYPE INTEGER DEFAULT 3 (The genes need to be expressed in at least this many cells.)
# PARAMETER OPTIONAL groupident: "Sample or group name" TYPE STRING DEFAULT empty (Type the group or sample name or identifier here. For example CTRL, STIM, TREAT. Do not use underscore _ in the names! Fill this field if you are combining samples later.)
# RUNTIME R-3.6.1


# 2017-06-06 ML
# 2017-07-05 ML split into separate tool
# 2018-01-11 ML update Seurat version to 2.2.0
# 2018-04-24 ML + AMS improve the input tar handling
# 2019-05-22 ML update Seurat version to 3.0
# 2019-09-24 ML correct for cases where there are NaNs in percent.mt
# 2019-09-30 EK add spport for lower case mitochondrial gene names

# Parameter removed from new R-version: "This functionality has been removed to simplify the initialization process/assumptions. 
# If you would still like to impose this threshold for your particular dataset, simply filter the input expression matrix before calling this function."
# PARAMETER OPTIONAL mingenes: "Keep cells which express at least this many genes" TYPE INTEGER DEFAULT 200 (The cells need to have expressed at least this many genes.)


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
	
	# Open tar package. Make a folder called datadir, open the tar there so that each file 
	# will be on the root level (remove everything from the name until the last "/" with the --xform option)
	system("mkdir datadir; cd datadir; tar xf ../files.tar --xform='s#^.+/##x' 2>> log.txt")	
	
	# Load the data
	dat <- Read10X("datadir/")	
	
}else{
	stop(paste('CHIPSTER-NOTE: ', "You need to provide either a 10X directory as a Tar package OR a DropSeq DGE as a tsv table. Please check your input file."))
}

# Initialize the Seurat object
seurat_obj <- CreateSeuratObject(counts = dat, min.cells = mincells, project = project.name)
# min.features = 200 => this is done in the next tool.
# v3: raw.data -> counts

# For sample detection later on
if (groupident != "empty") {
	seurat_obj@meta.data$stim <- groupident
}


# QC
# % of mito genes (note: they are named either "MT-CO1" or "mt-Co1", have to check both)
# NOTE: The pattern provided works for human and mouse gene names. You may need to adjust depending on your system of interest
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-|^mt-")


# pdf plots
pdf(file="QCplots.pdf", , width=13, height=7) 

# Violinplot
if (sum(is.na(seurat_obj@meta.data$percent.mt)) < 1) {
	VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
} else {
	VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
}

# FeatureScatter (v3)
plot1 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

# Number of cells:
textplot(paste("\v \v Number of \n \v \v cells: \n \v \v", length(colnames(x=seurat_obj))), halign="center", valign="center", cex=2)

dev.off() # close the pdf

# Save the Robj for the next tool
save(seurat_obj, file="seurat_obj.Robj")

## EOF

