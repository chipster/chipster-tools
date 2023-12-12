# TOOL single-cell-seurat-merge-samples.R: "Seurat v5 -Merge samples" (Merge multiple samples for joined analysis.)
# INPUT samples{...}.Robj: "Samples to combine" TYPE GENERIC
# OUTPUT seurat_obj_combined.Robj
# RUNTIME R-4.3.1-single-cell
# TOOLS_BIN ""

# 2023-11-29 IH

library(Seurat)
options(Seurat.object.assay.version = "v5")

# Load the R-Seurat-objects (called seurat_obj -that's why we need to rename them here)
input.names <- read.table("chipster-inputs.tsv", header = FALSE, sep = "\t")
for (i in 1:nrow(input.names)) {
    load(input.names[i, 1])
    name.of.obj <- paste("seurat_obj_", i, sep = "")
    assign(name.of.obj, seurat_obj)
}

seurat.objects.list <- as.list(mget(objects(pattern = "seurat_obj_")))
# merge first seurat object in the list with the other seurat object(s)
seurat_obj <- merge(x = seurat.objects.list[[1]], y = seurat.objects.list[-1])
seurat_obj

#moved from filter and regress tool
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj)

seurat_obj
# Save the Robj for the next tool
save(seurat_obj, file = "seurat_obj_combined.Robj")

## EOF
