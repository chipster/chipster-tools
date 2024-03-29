# TOOL single-cell-seurat-subset-v5.R: "Seurat v5 -Subset Seurat objects based on gene expression" (Subset cells in a Seurat object based on the expression level of a gene or feature. Gene name and expression threshold are given as parameters.)
# INPUT seurat_obj.Robj: "Seurat object" TYPE GENERIC
# OUTPUT seurat_obj_subset.Robj
# OUTPUT OPTIONAL log.txt
# OUTPUT OPTIONAL subset_plot.pdf
# PARAMETER OPTIONAL gene: "Gene" TYPE STRING DEFAULT "MS4A1" (Gene or feature name for subsetting. Set below the threshold for expression value.)
# PARAMETER OPTIONAL threshold: "Expression level threshold" TYPE DECIMAL DEFAULT 1 (Subset cells with higher than this expression in the gene selected above.)
# RUNTIME R-4.3.2-single-cell
# TOOLS_BIN ""

# 25.03.2020 ML
# 2021-10-04 ML Update to Seurat v4
# 2023-12-15 IH Update to Seurat v5

library(Seurat)
library(gplots)
options(Seurat.object.assay.version = "v5")

# Load the R-Seurat-object (called seurat_obj)
load("seurat_obj.Robj")

# Subset on the expression level of a gene/feature

# 1 sample:
if (exists("seurat_obj")) {
    # data.subset <- subset(x = seurat_obj, subset = gene > threshold)
    select.gene <- GetAssayData(object = seurat_obj, assay = "RNA", slot = "data")[gene, ]
    gene_ids <- names(which(select.gene > threshold))
    data.subset <- subset(seurat_obj, cells = gene_ids)
    # Save the subsetted Seurat Robj (rename to keep the object name for use in other tools)::
    seurat_obj <- data.subset
    save(seurat_obj, file = "seurat_obj_subset.Robj")
}

# 2 samples:
if (exists("data.combined")) {
    # data.subset <- subset(x = data.combined, subset = gene > threshold)
    select.gene <- GetAssayData(object = data.combined, assay = "RNA", slot = "data")[gene, ]
    gene_ids <- names(which(select.gene > threshold))
    data.subset <- subset(data.combined, cells = gene_ids)
    # Save the subsetted Seurat Robj (rename to keep the object name for use in other tools):
    data.combined <- data.subset
    save(data.combined, file = "seurat_obj_subset.Robj")
}

# Number of cells left, pdf plot:
pdf(file = "subset_plot.pdf", , width = 13, height = 7)
textplot(paste("\v \v Number of \n \v \v cells: \n \v \v", length(colnames(x = data.subset))), halign = "center", valign = "center", cex = 2)
# FeaturePlot(data.subset, gene) # would only work for data after clustering step
dev.off() # close the pdf

## EOF
