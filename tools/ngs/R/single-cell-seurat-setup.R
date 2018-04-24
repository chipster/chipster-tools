# TOOL single-cell-seurat-setup.R: "BETA Seurat -Setup and QC" (Setup the Seurat object, quality control, filter and regress the cells, determine statistically significant principal components. As an input, give a .tar package of a folder which contains the 10X output files OR a DGE matrix for DropSeq data. Please check that your input is assigned correctly under the parameters!)
# INPUT OPTIONAL files.tar: "tar package of 10X output files" TYPE GENERIC
# INPUT OPTIONAL dropseq.tsv: "DGE table from DropSeq" TYPE GENERIC
# OUTPUT OPTIONAL QCplots.pdf 
# OUTPUT OPTIONAL PCAplots.pdf
# OUTPUT OPTIONAL PCAgenes.txt
# OUTPUT OPTIONAL seurat_obj.Robj
# PARAMETER OPTIONAL project.name: "Project name for plotting" TYPE STRING DEFAULT Project_name (You can give your project a name. The name will appear on the plots.)
# PARAMETER OPTIONAL mincells: "Keep genes which are expressed in at least this many cells" TYPE INTEGER DEFAULT 3 (The genes need to be expressed in at least this many cells.)
# PARAMETER OPTIONAL mingenes: "Keep cells which express at least this many genes" TYPE INTEGER DEFAULT 200 (The cells need to have expressed at least this many genes.)
# RUNTIME R-3.4.3


# 2017-06-06 ML
# 2017-07-05 ML split into separate tool
# 2018-01-11 ML update Seurat version to 2.2.0
# 2018-04-24 ML + AMS improve the input tar handling

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


# Save the Robj for the next tool
save(seurat_obj, file="seurat_obj.Robj")

## EOF

