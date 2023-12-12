# TOOL single-cell-filter-v5.R: "Seurat v5 -Filter cells" (This tool filters out dead cells, empties and doublets. It then normalizes gene expression values and detects highly variable genes across the cells. Finally, it scales the data and regresses out unwanted variation based on the number of UMIs and mitochondrial transcript percentage. You can also choose to regress out variation due to cell cycle heterogeneity.)
# INPUT OPTIONAL seurat_obj.Robj: "Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL seurat_obj_filter.Robj
# PARAMETER OPTIONAL mingenes: "Filter out cells which have less than this many genes expressed" TYPE INTEGER DEFAULT 200 (Filter out empties. The cells to be kept must express at least this number of genes.)
# PARAMETER OPTIONAL genecountcutoff: "Filter out cells which have more than this many genes expressed" TYPE INTEGER DEFAULT 2500 (Filter out multiplets. The cells to be kept must express less than this number of genes.)
# PARAMETER OPTIONAL mitocutoff: "Filter out cells which have higher mitochondrial transcript percentage" TYPE DECIMAL FROM 0 TO 100 DEFAULT 5 (Filter out dead cells. The cells to be kept must have lower percentage of mitochondrial transcripts than this if needed in your data.)
# PARAMETER OPTIONAL minribo: "Filter out cells which have lower ribosomal transcript percentage" TYPE DECIMAL FROM 0 TO 100 DEFAULT 0 (Filter out cells that have lower ribosomal transcript percentage.)
# RUNTIME R-4.3.1-single-cell
# TOOLS_BIN ""


library(Seurat)
options(Seurat.object.assay.version = "v5")

# Load the R-Seurat-object (called seurat_obj)
load("seurat_obj.Robj")

# Subset: remove potential empties, multiplets and broken cells based on parameters.
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > mingenes & nFeature_RNA < genecountcutoff & percent.mt < mitocutoff & percent.rb >= minribo)

# Save the Robj for the next tool
save(seurat_obj, file = "seurat_obj_filter.Robj")

## EOF
