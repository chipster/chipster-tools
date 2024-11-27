# TOOL single-cell-seurat-separate-samples-from-sobject.R: "Seurat v5 BETA -Separate samples from Seurat object" (With this tool you can separate samples from a large Seurat object.)
# INPUT OPTIONAL seurat_obj.Robj: "Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL QCplots_{...}.pdf
# OUTPUT OPTIONAL seurat_obj_{...}.Robj
# RUNTIME R-4.3.2-single-cell
# SLOTS 5
# TOOLS_BIN ""

# 2024-11-25 ML

# RUNTIME R-4.2.3-seurat5
# RUNTIME R-4.2.3-single-cell

library(Seurat)
library(dplyr)
library(Matrix)
library(gplots)

load("seurat_obj.Robj")

# For testing:
# seurat_obj <- subset_of_samples

object_to_split <- seurat_obj

samples <- levels(object_to_split@meta.data$orig.ident)


length(samples)
how_many_samples <- length(samples)
# for testing
#  for (i in 1:3) {
for (i in 1:how_many_samples) {
    a_sample_name <- samples[i]
    a_sample <- subset(x = object_to_split, subset = orig.ident == a_sample_name)
    # print(a_sample_name)

    # Re-do the QC plots from the setup tool for each sample:
    pdf_name <- paste("QCplots_", a_sample_name, ".pdf", sep = "")
    pdf(file = pdf_name, , width = 13, height = 7)
    # Violinplot
    if ((sum(is.na(a_sample@meta.data$percent.mt)) < 1) && (sum(is.na(a_sample@meta.data$percent.rb)) < 1)) {
        print(VlnPlot(a_sample, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 4))
    } else if ((sum(is.na(a_sample@meta.data$percent.mt)) < 1)) {
        print(VlnPlot(a_sample, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
    } else if ((sum(is.na(a_sample@meta.data$percent.rb)) < 1)) {
        print(VlnPlot(a_sample, features = c("nFeature_RNA", "nCount_RNA", "percent.rb"), ncol = 3))
    } else {
        print(VlnPlot(a_sample, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2))
    }
    # Number of cells:
    textplot(paste("\v \v Number of \n \v \v cells: \n \v \v in sample ", a_sample_name, ": ", length(colnames(x = a_sample))), halign = "center", valign = "center", cex = 2)
    dev.off() # close the pdf

    # Save the subsetted Robj
    # rename as seurat_obj for later use in other tools:
    seurat_obj <- a_sample
    name.of.obj <- paste("seurat_obj_", a_sample_name, ".Robj", sep = "")
    save(seurat_obj, file = name.of.obj)
}

# EOF
