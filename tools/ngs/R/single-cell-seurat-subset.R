# TOOL single-cell-seurat-subset.R: "Seurat v3 BETA -Subset Seurat objects based on gene expression" (Subset cells in a Seurat object based on the expression level of a gene or feature. Gene name and expression threshold are given as parameters.)
# INPUT seurat_obj.Robj: "Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL log.txt
# OUTPUT OPTIONAL subset_plot.pdf
# OUTPUT seurat_obj_subset.Robj
# PARAMETER OPTIONAL gene: "Gene" TYPE STRING DEFAULT "MS4A1" (Gene or feature name for subsetting. Set below the threshold for expression value.)
# PARAMETER OPTIONAL threshold: "Expression level threshold" TYPE DECIMAL DEFAULT 1 (Subset cells with higher than this expression in the gene selected above.)
# RUNTIME R-3.6.1-single-cell

# 25.03.2020 ML

library(Seurat)
library(gplots)

# Load the R-Seurat-object (called seurat_obj)
load("seurat_obj.Robj")

# Subset on the expression level of a gene/feature

# 1 sample:
if (exists("seurat_obj")){
    # data.subset <- subset(x = seurat_obj, subset = gene > threshold)
    select.gene = GetAssayData(object = seurat_obj, assay = "RNA", slot = "data")[gene,]
    gene_ids = names(which(select.gene > threshold))
    data.subset = subset(seurat_obj, cells=gene_ids)
     # Save the subsetted Seurat Robj (rename to keep the object name for use in other tools)::
    seurat_obj <- data.subset
    save(seurat_obj, file="seurat_obj_subset.Robj")
}

# 2 samples:
if (exists("data.combined")){
    # data.subset <- subset(x = data.combined, subset = gene > threshold)
    select.gene = GetAssayData(object = data.combined, assay = "RNA", slot = "data")[gene,]
    gene_ids = names(which(select.gene > threshold))
    data.subset = subset(data.combined, cells=gene_ids)
    # Save the subsetted Seurat Robj (rename to keep the object name for use in other tools):
    data.combined <- data.subset
    save(data.combined, file="seurat_obj_subset.Robj")
}

# Number of cells left, pdf plot:
pdf(file="subset_plot.pdf", , width=13, height=7) 
textplot(paste("\v \v Number of \n \v \v cells: \n \v \v", length(colnames(x=data.subset))), halign="center", valign="center", cex=2)
# FeaturePlot(data.subset, gene) # would only work for data after clustering step
dev.off() # close the pdf

## EOF