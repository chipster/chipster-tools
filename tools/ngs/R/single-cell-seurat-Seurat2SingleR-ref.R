# TOOL single-cell-seurat-Seurat2SingleR-ref.R: "Seurat v5 - Create custom SingleR reference with pre-annotated seurat object"
# INPUT seurat_obj.Robj: "Seurat object. Has to be pre-processed so that it contains UMAP information." TYPE GENERIC
# INPUT seurat_obj_clustering.Robj: "Reference Seurat object with pre-annotated cell types." TYPE GENERIC
# OUTPUT OPTIONAL Plots.pdf
# RUNTIME R-4.5.1-seurat5
# TOOLS_BIN ""





#Iivari's SingleR functions:

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






library("Seurat")
library("SingleR")
library("scater")


load("seurat_obj.Robj", verbose = TRUE)

seurat_ref_obj <- seurat_obj


# Figure out how to load a second seurat obj
#load("seurat_obj_ref.Robj")

load("seurat_obj_clustering.Robj", verbose = TRUE)





ref <- build_singler_reference(seurat_obj, label_col = "seurat_clusters", aggr_ref = FALSE)

print("refdone")

# Unannotated seurat obj
expr <- GetAssayData(seurat_ref_obj, layer = "data")

expr <- as.matrix(expr)


anno <- SingleR(test = expr, ref = ref, assay.type.test=1,
                labels = ref$label)


# ══════════════════════════════════════════════════════════════════════════════
#  build_singler_reference()
#  Cell-level reference with optional per-label downsampling
# ══════════════════════════════════════════════════════════════════════════════

#Iivari func No. 2

#ref <- build_singler_reference(seurat_ref_obj, label_col = "seurat_clusters")

anno2 <- SingleR(test = expr, ref = ref, assay.type.test=1,
                 labels = ref$label)





seurat_obj <- run_singler_annotation(query_seurat = seurat_obj, ref = ref, label_col = "label")




seurat_obj$seurat$celltype <- seurat_obj$pred$labels

sobj <- seurat_obj$seurat

sobj$singler_label

sobj <- RunPCA(sobj, dims = 1:10)
sobj <- RunUMAP(sobj, dims = 1:10)


pdf(file = "Plots.pdf")

print(DimPlot(sobj, group.by = "singler_label"))

dev.off()


# EOF