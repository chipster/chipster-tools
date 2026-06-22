# TOOL single-cell-seurat-Seurat2SingleR-ref.R: "Seurat v5 - Annotate seurat object's celltypes with pre-annotated seurat object using SingleR"
# INPUT seurat_ref_obj.Robj: "Reference Seurat object with pre-annotated cell types. Object must be renamed as seurat_ref_obj.Robj" TYPE GENERIC
# INPUT seurat_obj_unannotated.Robj: "Seurat object. Has to be pre-processed so that it contains at least UMAP information. Has to be renamed as seurat_obj_unannotated.Robj" TYPE GENERIC
# OUTPUT seurat_obj_annotated.Robj
# OUTPUT Built_SingleR_ref.Robj
# OUTPUT OPTIONAL Plots.pdf
# PARAMETER OPTIONAL aggregate_reference: "Aggregate cells into one “pseudo-bulk” sample per label (e.g., by averaging across log-expression values) and using that as the reference profile. If set to TRUE, faster to run but may not be as accurate." TYPE [FALSE, TRUE] DEFAULT FALSE
# RUNTIME R-4.5.1-seurat5
# TOOLS_BIN ""




# 2026-06-22 JV


#  Iivari's SingleR functions:

#  build_singler_reference()
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
  if (is.null(labels)) stop("Column '", label_col, "' not found in metadata.")

  seurat_obj <- seurat_obj[!duplicated(rownames(seurat_obj)), ]

  
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



#  run_singler_annotation()
#  Annotate a query Seurat object using a SingleR reference
run_singler_annotation <- function(
    query_seurat,
    ref,
    label_col      = "label",
    assay          = "RNA",
    layer          = "data",
    fine_tune      = TRUE,
    prune_score    = TRUE,
    BPPARAM        = BiocParallel::MulticoreParam(workers = 4),
    add_to_seurat  = TRUE,
    result_prefix  = "singler"
) {
  
  mat <- LayerData(query_seurat, assay = assay, layer = layer)
  
  pred <- SingleR(
    test      = mat,
    ref       = ref,
    labels    = ref[[label_col]],
    fine.tune = fine_tune,
    prune     = prune_score,
    BPPARAM   = BPPARAM
  )
  
  if (add_to_seurat) {
    
    meta_df <- data.frame(
      row.names = rownames(pred),
      label     = pred$pruned.labels,
      score     = apply(pred$scores, 1, max, na.rm = TRUE),
      delta     = pred$delta.next
    )
    
    # Rename columns with prefix cleanly — no data.table := involved
    colnames(meta_df) <- paste0(result_prefix, "_", colnames(meta_df))
    
    query_seurat <- AddMetaData(query_seurat, metadata = meta_df)
    
    return(list(seurat = query_seurat, pred = pred))
  }
  
  pred
}



aggregate_reference <- as.logical(aggregate_reference)


library("Seurat")
library("SingleR")
library("scater")


load("seurat_ref_obj.Robj", verbose = TRUE)

seurat_ref_obj <- seurat_obj
rm(seurat_obj)

print(unique(colnames(seurat_ref_obj@meta.data)))

# Chipster pipeline puts celltypes as idents(??), needs to be added as a column for this function
seurat_ref_obj$celltype <- Idents(seurat_ref_obj)

head(seurat_ref_obj$celltype)
print("Up is the reference seurat object")



load("seurat_obj_unannotated.Robj", verbose = TRUE)

print(unique(colnames(seurat_obj@meta.data)))
print("Now this one is the unannotated sobj")


print("Starting to build the reference")


ref <- build_singler_reference(seurat_ref_obj, label_col = "celltype", aggr_ref = aggregate_reference)

print("Reference building succesful")

# Annotate the unannotated expression matrix by using the just created reference


print("Running singleR annotation on the seurat object")
seurat_obj <- run_singler_annotation(query_seurat = seurat_obj, ref = ref, label_col = "label")


seurat_obj$seurat$celltype <- seurat_obj$pred$labels

predictions <- seurat_obj$pred


sobj <- seurat_obj$seurat
sobj <- SetIdent(object = sobj, value = predictions$labels)

# New PCA and UMAP? Maybe needed, maybe not, check SingleR docs

#sobj <- RunPCA(sobj, dims = 1:10)
#sobj <- RunUMAP(sobj, dims = 1:10)

print("Annotation succesful")


pdf(file = "Plots.pdf")

print(DimPlot(sobj, group.by = "singler_label"))
print(plotScoreHeatmap(predictions))
print(plotDeltaDistribution(predictions))

dev.off()

save(ref, file = "Built_SingleR_ref.Robj")
save(sobj, file = "seurat_obj_annotated.Robj")

# EOF