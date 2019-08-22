# TOOL single-cell-seurat-update-object.R: "Seurat v3 -Update Seurat objects to version 3" (Updates Seurat objects to new v3 structure, so old objects can be analyzed with the updated versions of the tools.) 
# INPUT OPTIONAL seurat_obj.Robj: "Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL seurat_obj_2.Robj
# RUNTIME R-3.6.1

# 2019-07-24 ML 

library(Seurat)

# Load the R-Seurat-object (called seurat_obj)
load("seurat_obj.Robj")

seurat_obj_2.Robj = UpdateSeuratObject(object = seurat_obj)


# Save the Robj for the next tool
save(seurat_obj, file="seurat_obj_2.Robj")

## EOF