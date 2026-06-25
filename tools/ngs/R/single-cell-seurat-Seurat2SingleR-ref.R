# TOOL single-cell-seurat-Seurat2SingleR-ref.R: "Seurat v5 - Build celltype reference from Seurat object" (With this tool you can make a custom reference SummarizedExperiment object, which can be used to annotate cells in seurat objects)
# INPUT seurat_ref_obj.Robj: "Reference Seurat object with pre-annotated cell types." TYPE GENERIC
# OUTPUT SummarizedExperiment_reference.Robj
# PARAMETER OPTIONAL aggregate_reference: "Aggregate cells into one “pseudo-bulk” sample per label (e.g., by averaging across log-expression values) and using that as the reference profile. If set to TRUE, faster to run but may not be as accurate." TYPE [FALSE, TRUE] DEFAULT FALSE
# RUNTIME R-4.5.1-seurat5
# TOOLS_BIN ""

aggregate_reference <- as.logical(aggregate_reference)

# 2026-06-22 JV


#  Iivari's SingleR functions:

#  Create a SummarizedExperiment reference from a labelled Seurat object
build_singler_reference <- function(
    seurat_obj,
    label_col      = "celltype",   # metadata column with cell type labels
    assay          = "RNA",
    layer          = "data",       # use normalised counts
    genes          = NULL,         # optional: restrict to these genes
    aggr_ref       = TRUE,         # aggregate to pseudo-bulk per label
    aggr_args      = list(power = 0.5)
) {
  
  mat <- LayerData(seurat_obj, assay = assay, layer = layer)
  
  labels <- seurat_obj[[label_col, drop = TRUE]]
  if (is.null(labels)) stop("CHIPSTER-NOTE: Column '", label_col, "' not found in metadata.")

  
  if (!is.null(genes)) {
    genes <- intersect(genes, rownames(mat))
    mat   <- mat[genes, ]
  }
  
  ref <- SummarizedExperiment::SummarizedExperiment(
    assays = list(logcounts = mat),
    colData = S4Vectors::DataFrame(label = labels)
  )
  
  if (aggr_ref) {
    ref <- aggregateReference(ref, ref$label, power = aggr_args$power)
  }
  
  ref
}



library("Seurat")
library("SingleR")
library("scater")

# Load Seurat object that will be used as a reference

load("seurat_ref_obj.Robj", verbose = TRUE)

seurat_ref_obj <- seurat_obj
rm(seurat_obj)

# print(unique(colnames(seurat_ref_obj@meta.data)))

# Chipster pipeline puts celltypes as idents(??), needs to be added as a column for this function

if (any(Idents(seurat_ref_obj) %in% c(0,1,2))) {
  stop("CHIPSTER-NOTE: No cell types found, only cluster numbers.")
}

seurat_ref_obj$celltype <- Idents(seurat_ref_obj)

#head(seurat_ref_obj$celltype)




print("Starting to build the reference")


SummarizedExperiment_reference <- build_singler_reference(seurat_ref_obj, label_col = "celltype", aggr_ref = aggregate_reference)

print("Reference building succesful")


save(SummarizedExperiment_reference, file = "SummarizedExperiment_reference.Robj")


# EOF