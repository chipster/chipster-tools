# TOOL single-cell-seurat-split-samples-from-matrix.R: "Seurat v4 - BETA Separate samples" (Setup the Seurat object, make quality control plots and filter out genes. There are 3 options for input files, please check that your input file is correctly assigned under the parameters. If you have 10X data, make a tar package containing the files genes.tsv, barcodes.tsv and matrix.mtx \(you can use the tool \"Utilities - Make a tar package\" for this\). Alternatively you can give a DGE matrix or a 10X CellRanger hdf5 file as input. If you are planning to combine samples later on, make sure you name them in this tool!)
# INPUT OPTIONAL matrix.tsv: "Multisample DGE table in tsv format" TYPE GENERIC
# OUTPUT OPTIONAL sample_{...}.tsv
# PARAMETER OPTIONAL sample.name.to.look.for: "Sample name to separate into own matrix" TYPE STRING DEFAULT My_sample_name (Give here the name of the sample as it is in column names of the full matrix)
# RUNTIME R-4.2.3-single-cell
# SLOTS 2
# TOOLS_BIN ""


# 2023-08-24 ML


source(file.path(chipster.common.path, "tool-utils.R"))
#version <- system(paste(bowtie.binary,"--version | head -1 | cut -d ' ' -f 3"),intern = TRUE)
#package.version("Seurat")
#version <- package.version("Seurat")
#documentVersion("Seurat",version)

# If input is matrix:
if (file.exists("matrix.tsv")) {
  dat <- read.table("matrix.tsv", header = TRUE, sep = "\t", row.names = 1)
} 

# Remove empty spaces from user typed name:
sample.name.to.look.for <- gsub(" ", "", sample.name.to.look.for)

# Select the rows matching the sample name:
dat_sample <- dat[, grep(pattern =sample.name.to.look.for, colnames(dat))]   

# 
name.for.output.file <- paste("sample_", sample.name.to.look.for, ".tsv", sep="")


# Write to table
write.table(as.matrix(dat_sample), file = name.for.output.file, sep = "\t", row.names = T, col.names = T, quote = F)
  


## EOF


