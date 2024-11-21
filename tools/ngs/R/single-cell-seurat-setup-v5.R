# TOOL single-cell-seurat-setup-v5.R: "Seurat v5 -Setup and QC" (Setup the Seurat object, make quality control plots and filter out genes. There are several options for input files. Check that your input file is correctly assigned under the parameters. If you have 10X filtered feature-barcode matrix files in MEX format, make a tar package containing files genes.tsv, barcodes.tsv and matrix.mtx \(you can use the tool \"Utilities - Make a tar package\"\). Alternatively you can give a DGE matrix is tsv format, a 10X filtered feature-barcode matrix in hdf5 format or a CellBender filtered feature-barcode matrix in hdf5 format. If you are planning to combine samples later on, make sure you name them in this tool!)
# INPUT OPTIONAL files.tar: "tar package of 10X filtered feature-barcode matrix files in MEX format" TYPE GENERIC
# INPUT OPTIONAL dropseq.tsv: "DGE table in tsv format" TYPE GENERIC
# INPUT OPTIONAL hdf5.h5: "10X or CellBender filtered feature-barcode matrix in  hdf5 format" TYPE GENERIC
# OUTPUT OPTIONAL setup_seurat_obj.Robj
# OUTPUT OPTIONAL QCplots.pdf
# OUTPUT OPTIONAL PCAplots.pdf
# OUTPUT OPTIONAL PCAgenes.txt
# PARAMETER project.name: "Project name for plotting" TYPE STRING DEFAULT My_project_name (You can give your project a name. The name will appear on the plots. Do not use underscore _ in the names!)
# PARAMETER sample_name: "Sample name" TYPE STRING DEFAULT control1 (Type the sample name or identifier here. For example control1, cancer3a. Do not use underscore _ in the names! Fill this field if you are combining samples later.)
# PARAMETER sample.group: "Sample group" TYPE STRING DEFAULT CTRL (Type the sample name or identifier here. For example CTRL, STIM, TREAT. Do not use underscore _ in the names! Fill this field if you are combining samples later.)
# PARAMETER OPTIONAL mincells: "Keep genes which are expressed in at least this many cells" TYPE INTEGER DEFAULT 3 (The genes need to be expressed in at least this many cells.)
# RUNTIME R-4.3.2-single-cell
# SLOTS 8
# TOOLS_BIN ""



# 2017-06-06 ML
# 2017-07-05 ML split into separate tool
# 2018-01-11 ML update Seurat version to 2.2.0
# 2018-04-24 ML + AMS improve the input tar handling
# 2019-05-22 ML update Seurat version to 3.0
# 2019-09-24 ML correct for cases where there are NaNs in percent.mt
# 2019-09-30 EK add spport for lower case mitochondrial gene names
# 2021-10-04 ML Update to Seurat v4
# 2022-04-01 ML Add HDF5 input file option
# 2023-02-01 ML Add 5 slots
# 2023-04-06 LG Remove 5 slots
# 2023-02-01 ML Return to the original 2 slots
# 2023-06-14 ML Allow 10X tar input files in gzipped format and with longer file names
# 2023-07-12 ML Mitogenes with Mt-
# 2023-08-28 IH add percent.rb to QC
# 2023-10-10 IH Update to Seurat v5
# 2023-02-01 ML Slots from 2 -> 5


# Parameter removed from new R-version: "This functionality has been removed to simplify the initialization process/assumptions.
# If you would still like to impose this threshold for your particular dataset, simply filter the input expression matrix before calling this function."
# PARAMETER OPTIONAL mingenes: "Keep cells which express at least this many genes" TYPE INTEGER DEFAULT 200 (The cells need to have expressed at least this many genes.)


library(Seurat)
library(dplyr)
library(Matrix)
library(gplots)
library(Biobase)
options(Seurat.object.assay.version = "v5")
source(file.path(chipster.common.lib.path, "tool-utils.R"))
# version <- system(paste(bowtie.binary,"--version | head -1 | cut -d ' ' -f 3"),intern = TRUE)
package.version("Seurat")
version <- package.version("Seurat")
documentVersion("Seurat", version)

# Replace empty spaces in sample and project name with underscore _ :
sample_name <- gsub(" ", "_", sample_name)
project.name <- gsub(" ", "_", project.name)


# If using DropSeq data:
if (file.exists("dropseq.tsv")) {
  dat <- read.table("dropseq.tsv", header = TRUE, sep = "\t", row.names = 1)
  # Change the dataframe to a matrix
  dat <- Matrix(as.matrix(dat), sparse = TRUE)
  # If using 10X data:
} else if (file.exists("files.tar")) {
  # Read the contents of the tar file into a list
  system("tar tf files.tar > tar.contents 2>> log.txt")
  file.list <- scan("tar.contents", what = "", sep = "\n")

  # Check that the input is a valid tar file
  if (length(file.list) == 0) {
    stop(paste("CHIPSTER-NOTE: ", "It seems your input file is not a valid Tar package. Please check your input file."))
  }

  # Open tar package. Make a folder called datadir, open the tar there so that each file
  # will be on the root level (remove everything from the name until the last "/" with the --xform option)
  system("mkdir datadir; cd datadir; tar xf ../files.tar --xform='s#^.+/##x' 2>> log.txt")

  # If the features.tsv, barcodes.tsv and matrix.mtx files are gzipped, open them:
  if (length(list.files(path = "datadir/", pattern = "*matrix.mtx.gz")) >= 1) {
    system("cd datadir; gzip -d *matrix.mtx.gz")
  }
  if (length(list.files(path = "datadir/", pattern = "*features.tsv.gz")) >= 1) {
    system("cd datadir; gzip -d *features.tsv.gz")
  }
  if (length(list.files(path = "datadir/", pattern = "*genes.tsv.gz")) >= 1) {
    system("cd datadir; gzip -d *genes.tsv.gz")
  }
  if (length(list.files(path = "datadir/", pattern = "*barcodes.tsv.gz")) >= 1) {
    system("cd datadir; gzip -d *barcodes.tsv.gz")
  }

  # If features.tsv, barcodes.tsv and matrix.mtx files have something extra in the name, remove the extra:
  system("cd datadir; mv *features.tsv features.tsv")
  system("cd datadir; mv *genes.tsv genes.tsv")
  system("cd datadir; mv *barcodes.tsv barcodes.tsv")
  system("cd datadir; mv *matrix.mtx matrix.mtx")

  # rename features.tsv as genes.tsv (Read10X SHOULD work with features.tsv as well, but it doesn't, yet at least)
  if (file.exists("datadir/features.tsv")) {
    system("mv datadir/features.tsv datadir/genes.tsv")
  }

  # Load the data
  dat <- Read10X("datadir/")

  # If using HDF5 data:
} else if (file.exists("hdf5.h5")) {
  library(hdf5r)
  dat <- Read10X_h5(filename = "hdf5.h5", use.names = T)
} else {
  stop(paste("CHIPSTER-NOTE: ", "You need to provide either a 10X directory as a Tar package OR a DropSeq DGE as a tsv table OR a 10X h5 file. Please check your input file."))
}

# Initialize the Seurat object

# In case of multiomics data, "dat" object has a list of matrices.
# We need to select the Gene Expression one from there.
if (class(dat) == "list") {
  seurat_obj <- CreateSeuratObject(counts = dat$`Gene Expression`, min.cells = mincells, project = project.name)
} else {
  seurat_obj <- CreateSeuratObject(counts = dat, min.cells = mincells, project = project.name)
  # min.features = 200 => this is done in the next tool.
}

# For sample detection later on
seurat_obj@meta.data$stim <- sample_name
seurat_obj$type <- sample.group

# QC
# % of mito genes (note: they are named either "MT-CO1" or "mt-Co1", have to check both)
# NOTE: The pattern provided works for human and mouse gene names. You may need to adjust depending on your system of interest
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-|^mt-|^Mt-")
seurat_obj[["percent.rb"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RPS|^RPL|^rps|^rpl|^Rps|^Rpl")


# pdf plots
pdf(file = "QCplots.pdf", , width = 13, height = 7)

# Violinplot
if ((sum(is.na(seurat_obj@meta.data$percent.mt)) < 1) && (sum(is.na(seurat_obj@meta.data$percent.rb)) < 1)) {
  VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 4)
} else if ((sum(is.na(seurat_obj@meta.data$percent.mt)) < 1)) {
  VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
} else if ((sum(is.na(seurat_obj@meta.data$percent.rb)) < 1)) {
  VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.rb"), ncol = 3)
} else {
  VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
}

# FeatureScatter (v3)
plot1 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

# Number of cells:
textplot(paste("\v \v Number of \n \v \v cells: \n \v \v", length(colnames(x = seurat_obj))), halign = "center", valign = "center", cex = 2)

dev.off() # close the pdf

# Save the Robj for the next tool
save(seurat_obj, file = "setup_seurat_obj.Robj")

## EOF
