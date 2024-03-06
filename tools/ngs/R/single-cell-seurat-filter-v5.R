# TOOL single-cell-seurat-filter-v5.R: "Seurat v5 -Filter cells" (This tool filters out dead cells, empties and doublets.)
# INPUT OPTIONAL seurat_obj.Robj: "Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL seurat_obj_filter.Robj
# OUTPUT OPTIONAL QCplots.pdf
# PARAMETER OPTIONAL mingenes: "Filter out cells which have less than this many genes expressed" TYPE INTEGER DEFAULT 200 (Filter out empties. The cells to be kept must express at least this number of genes.)
# PARAMETER OPTIONAL genecountcutoff: "Filter out cells which have more than this many genes expressed" TYPE INTEGER DEFAULT 2500 (Filter out multiplets. The cells to be kept must express less than this number of genes.)
# PARAMETER OPTIONAL mitocutoff: "Filter out cells which have higher mitochondrial transcript percentage" TYPE DECIMAL FROM 0 TO 100 DEFAULT 5 (Filter out dead cells. The cells to be kept must have lower percentage of mitochondrial transcripts than this if needed in your data.)
# PARAMETER OPTIONAL minribo: "Filter out cells which have lower ribosomal transcript percentage" TYPE DECIMAL FROM 0 TO 100 DEFAULT 0 (Filter out cells that have lower ribosomal transcript percentage.)
# RUNTIME R-4.3.2-single-cell
# TOOLS_BIN ""


library(Seurat)
library(dplyr)
library(Matrix)
library(gplots)
options(Seurat.object.assay.version = "v5")

# Load the R-Seurat-object (called seurat_obj)
load("seurat_obj.Robj")

# Subset: remove potential empties, multiplets and broken cells based on parameters.
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA >= mingenes & nFeature_RNA <= genecountcutoff & percent.mt <= mitocutoff & percent.rb >= minribo)


# Re-do the QC plots from the setup tool, for comparison:
# pdf plots
pdf(file = "QCplots.pdf", , width = 13, height = 7)

# Switch back to orig.ident from cell cycle scores:
# Idents(object = seurat_obj) <- "orig.ident"

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

# Number of cells:
textplot(paste("\v \v Number of \n \v \v cells: \n \v \v", length(colnames(x = seurat_obj))), halign = "center", valign = "center", cex = 2)

dev.off() # close the pdf

# Save the Robj for the next tool
save(seurat_obj, file = "seurat_obj_filter.Robj")

## EOF
