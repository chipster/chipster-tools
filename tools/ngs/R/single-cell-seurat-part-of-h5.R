# TOOL single-cell-seurat-part-of-h5.R: "Seurat v4 - BETA Extract only gene expression data from a multiomics h5 object" (In case you have a multiomics h5 object, you can use this tool to only get the gene expression matrix in tsv format. This matrix can be used as an input for the setup tool.)
# INPUT OPTIONAL hdf5.h5: "10X CellRanger hdf5 input file" TYPE GENERIC
# OUTPUT OPTIONAL sample_{...}.tsv
# PARAMETER OPTIONAL sample.name.for.file: "Output file name" TYPE STRING DEFAULT My_sample_name (Give here a name for the output file.)
# RUNTIME R-4.2.3-single-cell
# SLOTS 5
# TOOLS_BIN ""

# 2023-10-17 ML

library(hdf5r)
library(Seurat)

# source(file.path(chipster.common.path, "tool-utils.R"))
# version <- system(paste(bowtie.binary,"--version | head -1 | cut -d ' ' -f 3"),intern = TRUE)
# package.version("Seurat")
# version <- package.version("Seurat")
# documentVersion("Seurat",version)

# Load the data:
if (file.exists("hdf5.h5")) {
  library(hdf5r)
  dat <- Read10X_h5(filename = "hdf5.h5", use.names = T)
} else {
  stop(paste("CHIPSTER-NOTE: ", "You need to provide a multiomics h5 file as input. Please check your input file."))
}

# Choose only the gene expression part:
data_gene_expr <- dat[["Gene Expression"]]


# Write to table:
name.for.output.file <- paste("sample_", sample.name.for.file, ".tsv", sep = "")
write.table(as.matrix(data_gene_expr), file = name.for.output.file, sep = "\t", row.names = T, col.names = T, quote = F)


## EOF