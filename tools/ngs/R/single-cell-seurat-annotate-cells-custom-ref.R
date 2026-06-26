# TOOL single-cell-seurat-annotate-cells-custom-ref.R: "Seurat v5 - Annotate cells with own reference" (Annotate cells on a Seurat object by using your custom SummarizedExperiment object as the reference.)
# INPUT SummarizedExperiment_reference.Robj: "Reference object" TYPE GENERIC (A SummarizedExperiment object.)
# INPUT seurat_obj_unannotated.Robj: "Seurat object that will get annotated" TYPE GENERIC (Has to be pre-processed so that it contains at least UMAP information.)
# OUTPUT seurat_obj_annotated.Robj
# OUTPUT SingleR_custom_annotation_Plots.pdf
# PARAMETER OPTIONAL prune: "Pruning" TYPE [FALSE: "no", TRUE: "yes"] DEFAULT TRUE (If yes, removes weak cell types and will be set as NA.) 
# PARAMETER OPTIONAL fine.tune: "Fine tuning" TYPE [FALSE: "no", TRUE: "yes"] DEFAULT TRUE (If yes, improves ranking accuracy of the best label.) 
# PARAMETER OPTIONAL label.size: "Label size in the output plots" TYPE DECIMAL DEFAULT 4 (Label size for cluster numbers or cell type names on top of UMAP. If you don't want any labels, set this to 0.)
# PARAMETER OPTIONAL width: "Width of the output plots" TYPE INTEGER DEFAULT 10 (Width of the output plots in inches.)
# PARAMETER OPTIONAL height: "Height of the output plots" TYPE INTEGER DEFAULT 10 (Height of the output plots in inches.)
# RUNTIME R-4.5.1-seurat5
# TOOLS_BIN ""

chipster.threads.max <- as.numeric(chipster.threads.max)
prune <- as.logical(prune)
fine.tune <- as.logical(fine.tune)



# Function for SingleR annotation by Iivari Kleino


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



load("SummarizedExperiment_reference.Robj")


SummarizedExperiment_reference

# The actual R variable has to be named as SummarizedExperiment_refernce (Chipster does this with Build celltype ref form seurat object)
# This if exists is basically to check whether user actually inputted the correct file (Note that currently if they input a SummarizedExperiment object that is not named as stated before, this error will pop out)

if (!exists("SummarizedExperiment_reference")) {
  stop("CHIPSTER-NOTE: Wrong input file. Try swapping input files ")
}

ref <- SummarizedExperiment_reference
rm(SummarizedExperiment_reference)



load("seurat_obj_unannotated.Robj", verbose = TRUE)



print("Running SingleR annotation on the unannotated Seurat object")
seurat_obj <- run_singler_annotation(query_seurat = seurat_obj, ref = ref, label_col = "label")


# Add labels into the seurat object
seurat_obj$seurat$celltype <- seurat_obj$pred$labels

predictions <- seurat_obj$pred


seurat_obj <- seurat_obj$seurat
seurat_obj <- SetIdent(object = seurat_obj, value = predictions$labels)

# Assign one cell type per cluster 

seurat_table <- table(seurat_obj$seurat_clusters, seurat_obj$celltype)

type <- apply(seurat_table, 1, function(x) names(which.max(x)))


# Check that cluster numbers on the seurat object and on the table actually match


match <- all(seurat_obj$seurat_clusters == names(type[as.character(seurat_obj$seurat_clusters)]))

if (!match) {
  "CHIPSER-NOTE: Cluster number on the seurat object does not match with the seurat_table. Wrong cell types would be annotated. Stopping process..."
} else {
  print("Cluster numbers match, assigning a cell type per cluster")
}

# head(type, 20)
seurat_obj$cluster_celltype <- as.vector(type[as.character(seurat_obj$seurat_clusters)])

print("Annotation succesful")

if (length(predictions$pruned.labels) > 0) {
print("Pruned, saving QC plots")

pdf(file = "SingleR_custom_annotation_Plots.pdf", width = width, height = height)

p0 <- DimPlot(seurat_obj, group.by = "singler_label", label = T, label.size = label.size)
print(p0)

p1 <- plotScoreHeatmap(predictions)
print(p1)

p2 <- plotDeltaDistribution(predictions)
print(p2)

p3 <- DimPlot(seurat_obj, group.by = "cluster_celltype", label = T, label.size = label.size)

p4 <-DimPlot(seurat_obj, group.by = "seurat_clusters", label = T, label.size = label.size)

print(p3+p4)

dev.off()
save(seurat_obj, file = "seurat_obj_annotated.Robj")

} else {

print("Not pruned, no QC plots available.")
pdf(file = "SingleR_custom_annotation_Plots.pdf", width = width, height = height)
p0 <- DimPlot(seurat_obj, group.by = "singler_label", label = T, label.size = label.size)

print(p0)

p3 <- DimPlot(seurat_obj, group.by = "cluster_celltype", label = T, label.size = label.size)

p4 <- DimPlot(seurat_obj, group.by = "seurat_clusters", label = T, label.size = label.size)

print(p3+p4)

dev.off()

save(seurat_obj, file = "seurat_obj_annotated.Robj")
}



# EOF