# TOOL single-cell-seurat-Seurat2SingleR-ref.R: "Seurat v5 - Build a celltype reference from a pre-annotated Seurat and transfer labels to unannotated Seurat"
# INPUT seurat_ref_obj.Robj: "Reference Seurat object with pre-annotated cell types." TYPE GENERIC
# INPUT OPTIONAL seurat_obj_unannotated.Robj: "Seurat object. Has to be pre-processed so that it contains at least UMAP information." TYPE GENERIC
# OUTPUT OPTIONAL seurat_obj_annotated.Robj
# OUTPUT Built_SingleR_ref.Robj
# OUTPUT OPTIONAL Plots.pdf
# PARAMETER OPTIONAL aggregate_reference: "Aggregate cells into one “pseudo-bulk” sample per label (e.g., by averaging across log-expression values) and using that as the reference profile. If set to TRUE, faster to run but may not be as accurate." TYPE [FALSE, TRUE] DEFAULT FALSE
# PARAMETER OPTIONAL prune: "Pruning" TYPE [FALSE, TRUE] DEFAULT TRUE (If set to TRUE, removes weak cell types and will be set as NA.) 
# PARAMETER OPTIONAL fine.tune: "fine tuning" TYPE [FALSE, TRUE] DEFAULT TRUE (If set to TRUE, improves ranking accuracy of the best label.) 
# RUNTIME R-4.5.1-seurat5
# TOOLS_BIN ""

fine.tune <- as.logical(fine.tune)
prune <- as.logical(prune)
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
  if (is.null(labels)) stop("Column '", label_col, "' not found in metadata.")

  
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


#  Annotate a query Seurat object using a SingleR reference

run_singler_annotation <- function(
    query_seurat,
    ref,
    label_col      = "label",
    assay          = "RNA",
    layer          = "data",
    fine_tune      = fine.tune,
    prune_score    = prune,
    BPPARAM        = BiocParallel::MulticoreParam(workers = chipster.threads.max),
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
      label     = if (prune) pred$pruned.labels else pred$labels,
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

# Load Seurat object that will be used as a reference

load("seurat_ref_obj.Robj", verbose = TRUE)

seurat_ref_obj <- seurat_obj
rm(seurat_obj)

# print(unique(colnames(seurat_ref_obj@meta.data)))

# Chipster pipeline puts celltypes as idents(??), needs to be added as a column for this function

if (any(Idents(seurat_ref_obj) %in% c(0,1,2))) {
  stop("CHIPSTER-NOTE: No cell types found, only cluster numbers. \n Try swapping input files")
}

seurat_ref_obj$celltype <- Idents(seurat_ref_obj)

#head(seurat_ref_obj$celltype)




print("Starting to build the reference")


ref <- build_singler_reference(seurat_ref_obj, label_col = "celltype", aggr_ref = aggregate_reference)

print("Reference building succesful")


save(ref, file = "Built_SingleR_ref.Robj")


if (file.exists("seurat_obj_unannotated.Robj")) {

load("seurat_obj_unannotated.Robj", verbose = TRUE)

#print(unique(colnames(seurat_obj@meta.data)))


# Annotate the unannotated expression matrix by using the just created reference

print("Running SingleR annotation on the unannotated Seurat object")
seurat_obj <- run_singler_annotation(query_seurat = seurat_obj, ref = ref, label_col = "label")


# Add labels into the seurat object
seurat_obj$seurat$celltype <- seurat_obj$pred$labels

predictions <- seurat_obj$pred


sobj <- seurat_obj$seurat
sobj <- SetIdent(object = sobj, value = predictions$labels)

print("Annotation succesful")


pdf(file = "Plots.pdf", onefile = T)

p1 <- (DimPlot(sobj, group.by = "singler_label"))
print(p1)

print("If prune is set to TRUE, ScoreHeatmap & DeltaDistribution plots will be plotted")

if (prune) {

p2 <- (plotScoreHeatmap(predictions))
print(p2)

p3 <- (plotDeltaDistribution(predictions))
print(p3)

}

dev.off()

save(sobj, file = "seurat_obj_annotated.Robj")

}

# EOF