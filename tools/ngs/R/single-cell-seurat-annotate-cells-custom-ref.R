# TOOL single-cell-seurat-annotate-cells-custom-ref.R: "Seurat v5 - Annotate cells with custom SummarizedExperiment reference" (Annotate cells on a Seurat object by using your custom SummarizedExperiment object as the reference.)
# INPUT custom_singleR_ref.Robj: "Reference object" TYPE GENERIC
# INPUT seurat_obj_unannotated.Robj: "Seurat object. Has to be pre-processed so that it contains at least UMAP information." TYPE GENERIC
# OUTPUT seurat_obj_annotated.Robj
# OUTPUT Plots.pdf
# PARAMETER OPTIONAL prune: "Pruning" TYPE [FALSE, TRUE] DEFAULT TRUE (If set to TRUE, removes weak cell types and will be set as NA.) 
# PARAMETER OPTIONAL fine.tune: "fine tuning" TYPE [FALSE, TRUE] DEFAULT TRUE (If set to TRUE, improves ranking accuracy of the best label.) 
# RUNTIME R-4.5.1-seurat5
# TOOLS_BIN ""

chipster.threads.max <- as.numeric(chipster.threads.max)
prune <- as.logical(prune)
fine.tune <- as.logical(fine.tune)



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



load("custom_singleR_ref.Robj")

if (!exists("custom_singleR_ref")) {
stop("CHIPSTER-NOTE: Wrong input file. Try swapping input files ")
}

ref <- custom_singleR_ref
rm(custom_singleR_ref)



load("seurat_obj_unannotated.Robj", verbose = TRUE)



print("Running SingleR annotation on the unannotated Seurat object")
seurat_obj <- run_singler_annotation(query_seurat = seurat_obj, ref = ref, label_col = "label")


# Add labels into the seurat object
seurat_obj$seurat$celltype <- seurat_obj$pred$labels

predictions <- seurat_obj$pred


seurat_obj <- seurat_obj$seurat
seurat_obj <- SetIdent(object = seurat_obj, value = predictions$labels)

print("Annotation succesful")

if (length(predictions$pruned.labels) > 0) {
print("Pruned, saving QC plots")

pdf(file = "Plots.pdf")

p0 <- DimPlot(seurat_obj, group.by = "singler_label")
print(p0)

p1 <- plotScoreHeatmap(predictions)
print(p1)

p2 <- plotDeltaDistribution(predictions)
print(p2)

dev.off()
save(seurat_obj, file = "seurat_obj_annotated.Robj")

} else {

pdf(file = "Plots.pdf")
p0 <- DimPlot(seurat_obj, group.by = "singler_label")

print(p0)

dev.off()

save(seurat_obj, file = "seurat_obj_annotated.Robj")
}



# EOF