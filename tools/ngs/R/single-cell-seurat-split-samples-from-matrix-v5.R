# TOOL single-cell-seurat-split-samples-from-matrix-v5.R: "Seurat v5 - BETA Separate samples" (In case you have a large digital expression matrix containing multiple samples in one .tsv table, where the sample names are included in the column name together with the cell barcode, you can use this tool to separate the sample of a specific name into a separate .tsv table. This .tsv table can then be used as a input for the Setup tool.)
# INPUT OPTIONAL matrix.tsv: "Multisample DGE table in tsv format" TYPE GENERIC
# OUTPUT OPTIONAL sample_{...}.tsv
# PARAMETER OPTIONAL sample.name.to.look.for: "Sample name to separate into own matrix" TYPE STRING DEFAULT My_sample_name (Give here the name of the sample as it is in column names of the full matrix)
# RUNTIME R-4.3.2-single-cell
# SLOTS 2
# TOOLS_BIN ""


# 2023-08-24 ML
# 2023-12-15 IH


source(file.path(chipster.common.lib.path, "tool-utils.R"))
options(Seurat.object.assay.version = "v5")
# version <- system(paste(bowtie.binary,"--version | head -1 | cut -d ' ' -f 3"),intern = TRUE)
# package.version("Seurat")
# version <- package.version("Seurat")
# documentVersion("Seurat",version)

# If input is matrix:
if (file.exists("matrix.tsv")) {
  dat <- read.table("matrix.tsv", header = TRUE, sep = "\t", row.names = 1)
}

# Remove empty spaces from user typed name:
sample.name.to.look.for <- gsub(" ", "", sample.name.to.look.for)

# Select the rows matching the sample name:
dat_sample <- dat[, grep(pattern = sample.name.to.look.for, colnames(dat))]

#
name.for.output.file <- paste("sample_", sample.name.to.look.for, ".tsv", sep = "")


# Write to table
write.table(as.matrix(dat_sample), file = name.for.output.file, sep = "\t", row.names = T, col.names = T, quote = F)



## EOF
