# TOOL single-cell-seurat-annotate-cells-custom-ref.R: "Seurat v5 - Annotate cells with SingleR and own reference" (Annotate cells on a Seurat object by using your custom SummarizedExperiment object as the reference.)
# INPUT seurat_obj_unannotated.Robj: "Seurat object that will get annotated" TYPE GENERIC (Has to be pre-processed so that it contains at least UMAP information.)
# INPUT SummarizedExperiment_reference.Robj: "Reference object" TYPE GENERIC (A SummarizedExperiment object.)
# OUTPUT seurat_obj_custom_ref_annotated.Robj
# OUTPUT SingleR_custom_ref_annotation_Plots.pdf
# OUTPUT OPTIONAL cluster_celltype_table.tsv
# PARAMETER method: "Method to assign one cell type per cluster" TYPE [majority: "Majority", SingleR: "SingleR"] DEFAULT SingleR (SingleR computes an aggregated profile per cluster. Majority method checks which type of cells are found the most in a cluster and assigns that cell type for that cluster.)
# PARAMETER OPTIONAL prune: "Pruning" TYPE [FALSE: "no", TRUE: "yes"] DEFAULT TRUE (Remove weak cell types and set them as NA.) 
# PARAMETER OPTIONAL fine.tune: "Fine tuning" TYPE [FALSE: "no", TRUE: "yes"] DEFAULT TRUE (Improve ranking accuracy of the best label.) 
# PARAMETER OPTIONAL height: "Height of the output plots" TYPE INTEGER DEFAULT 10 (Height of the output plots in inches.)
# PARAMETER OPTIONAL width: "Width of the output plots" TYPE INTEGER DEFAULT 10 (Width of the output plots in inches.)
# PARAMETER OPTIONAL label.size: "Label size in the output plots" TYPE DECIMAL DEFAULT 4 (Label size for cluster numbers or cell type names on top of UMAP. If you don't want any labels, set this to 0.)
# RUNTIME R-4.5.1-seurat5
# SLOTS 2
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
    BPPARAM   = BPPARAM,
    de.method = "wilcox"
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

# 27-08-2026 JV Add majority and SingleR cluster annotation options

# Load libraries

library("Seurat")
library("SingleR")
library("scater")



load("SummarizedExperiment_reference.Robj")


# The actual ref object R variable has to be named as SummarizedExperiment_reference (Chipster does this with Build celltype ref form seurat object)

if (!exists("SummarizedExperiment_reference")) {
  stop("CHIPSTER-NOTE: Wrong input file. Try swapping input files ")
}

ref <- SummarizedExperiment_reference
rm(SummarizedExperiment_reference)



load("seurat_obj_unannotated.Robj")

if (exists("data.combined")) {
  seurat_obj <- data.combined
}

assay <- as.character(DefaultAssay(seurat_obj))


print("Running SingleR annotation on the unannotated Seurat object")
seurat_obj <- run_singler_annotation(query_seurat = seurat_obj, ref = ref, label_col = "label", assay = assay)


# Add labels into the seurat object
seurat_obj$seurat$celltype <- seurat_obj$pred$labels

predictions <- seurat_obj$pred


seurat_obj <- seurat_obj$seurat
seurat_obj <- SetIdent(object = seurat_obj, value = predictions$labels)



# Assign one cell type per cluster 

# Method == majority, cluster cell type assigned with the majority vote

if (method == "majority") {
seurat_table <- table(seurat_obj$seurat_clusters, seurat_obj$celltype)

type <- apply(seurat_table, 1, function(x) names(which.max(x)))


# Check that cluster numbers on the seurat object and on the table actually match
match <- all(seurat_obj$seurat_clusters == names(type[as.character(seurat_obj$seurat_clusters)]))

if (!match) {
  stop("CHIPSER-NOTE: Cluster number on the seurat object does not match with the seurat_table. Wrong cell types would be annotated. Stopping process...")
} else {
  print("Cluster numbers match, assigning a cell type per cluster")
}

seurat_obj$cluster_celltype <- as.vector(type[as.character(seurat_obj$seurat_clusters)])

print("Annotation succesful")

if (length(predictions$pruned.labels) > 0) {
print("Pruned, saving QC plots")

pdf(file = "SingleR_custom_ref_annotation_Plots.pdf", width = width, height = height)

p0 <- DimPlot(seurat_obj, group.by = "singler_label", label = T, label.size = label.size)
print(p0)

p1 <- plotScoreHeatmap(predictions)
print(p1)

p2 <- plotDeltaDistribution(predictions)
print(p2)

p3 <- DimPlot(seurat_obj, group.by = "cluster_celltype", label = T, label.size = label.size)

print(p3)

dev.off()
save(seurat_obj, file = "seurat_obj_custom_ref_annotated.Robj")
write.table(seurat_table, file = "cluster_celltype_table.tsv", sep = "\t", row.names = T, col.names = T, quote = F)

} else {

print("Not pruned, no QC plots available.")
pdf(file = "SingleR_custom_ref_annotation_Plots.pdf", width = width, height = height)
p0 <- DimPlot(seurat_obj, group.by = "singler_label", label = T, label.size = label.size)

print(p0)

p3 <- DimPlot(seurat_obj, group.by = "cluster_celltype", label = T, label.size = label.size)

print(p3)

dev.off()

save(seurat_obj, file = "seurat_obj_custom_ref_annotated.Robj")
write.table(seurat_table, file = "cluster_celltype_table.tsv", sep = "\t", row.names = T, col.names = T, quote = F)
 }
} 


# Method == singleR (SingleR function assings cell type for cluster)

if (method == "SingleR") {

Idents(seurat_obj) <- seurat_obj@meta.data$seurat_clusters

sce <- as.SingleCellExperiment(seurat_obj)

annotation.main.clusters <- SingleR(test = sce, assay.type.test = 1, ref = ref, labels = ref$label, clusters = Idents(seurat_obj), fine.tune = fine.tune, prune = prune)
new.cluster.ids <- annotation.main.clusters$labels
names(new.cluster.ids) <- levels(seurat_obj)
seurat_obj_annotated_clusters <- seurat_obj
seurat_obj_annotated_clusters <- RenameIdents(seurat_obj, new.cluster.ids)


print("Annotation succesful")

if (length(predictions$pruned.labels) > 0) {
print("Pruned, saving QC plots")

pdf(file = "SingleR_custom_ref_annotation_Plots.pdf", width = width, height = height)

p0 <- DimPlot(seurat_obj, group.by = "singler_label", label = T, label.size = label.size)
print(p0)

p1 <- plotScoreHeatmap(predictions)
print(p1)

p2 <- plotDeltaDistribution(predictions)
print(p2)

p3 <- DimPlot(seurat_obj_annotated_clusters, label = T, label.size = label.size) + ggtitle("Main annotations based on clusters")

print(p3)

dev.off()


seurat_obj$cluster_celltype <- as.character(Idents(seurat_obj_annotated_clusters))[colnames(seurat_obj)]
save(seurat_obj, file = "seurat_obj_custom_ref_annotated.Robj")

} else {

print("Not pruned, no QC plots available.")

pdf(file = "SingleR_custom_ref_annotation_Plots.pdf", width = width, height = height)
p0 <- DimPlot(seurat_obj, group.by = "singler_label", label = T, label.size = label.size)

print(p0)

p3 <- DimPlot(seurat_obj_annotated_clusters, label = T, label.size = label.size) + ggtitle("Main annotations based on clusters")

print(p3)

dev.off()

seurat_obj$cluster_celltype <- as.character(Idents(seurat_obj_annotated_clusters))[colnames(seurat_obj)]
save(seurat_obj, file = "seurat_obj_custom_ref_annotated.Robj")

 }
}
# EOF
