# TOOL spatial-transcriptomics-seurat-filter-v5-HD.R: "Seurat v5 HD -Filter spots" (This tool filters out spots based on mitochondrial, ribosomal and hemoglobin transcript percentage.)
# INPUT OPTIONAL seurat_spatial_setup.Robj: "Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL seurat_obj_filtered.Robj
# PARAMETER OPTIONAL mitocutoff: "Filter out spots which have higher mitochondrial transcript percentage" TYPE DECIMAL FROM 0 TO 100 DEFAULT 20 (Filter out spots that have higher percentage of mitochondrial transcripts than this.)
# PARAMETER OPTIONAL minribo: "Filter out spots which have lower ribosomal transcript percentage" TYPE DECIMAL FROM 0 TO 100 DEFAULT 0 (Filter out spots that have lower ribosomal transcript percentage than this.)
# PARAMETER OPTIONAL hbcutoff: "Filter out spots which have higher hemoglobin transcript percentage" TYPE DECIMAL FROM 0 TO 100 DEFAULT 20 (Filter out spots that have higher percentage of hemoglobin transcripts than this.)
# RUNTIME R-4.5.1-seurat5
# SLOTS 3
# TOOLS_BIN ""


# PARAMETER OPTIONAL genes: "Genes to filter" TYPE [NULL: Null, mt: Mitochondrial, Hb: Hemoglobin, both: Both] DEFAULT NULL (You can choose to remove genes based on the top expressed genes plot for example mitochondrial or hemoglobin genes.)

# 2026-02 ML

library(Seurat)
library(Biobase)

print(package.version("Seurat"))
documentVersion("Seurat", package.version("Seurat"))


# Load the R-Seurat-object (called seurat_obj)
load("seurat_spatial_setup.Robj")
print(Version(seurat_obj))

## PARAMETER OPTIONAL which_bin: "Filter based on this bin" TYPE [008um: 8um, 016um: 16um] DEFAULT 008um (Which bin is used for filtering according the thresholds set below.)
## Choose the bin (based on which_bin parameter)
# DefaultAssay(seurat_obj) <- paste("Spatial.", which_bin, sep="") # "Spatial.008um"
# nCount_bin <- paste("nCount_Spatial.", which_bin, sep="")
# nFeature_bin <- paste("nFeature_Spatial.", which_bin, sep="")

# Subset: remove potential empties, multiplets and broken cells based on parameters
seurat_obj <- subset(seurat_obj, subset = percent.mt <= mitocutoff & percent.hb <= hbcutoff & percent.rb >= minribo)

# Save the Robj for the next tool
save(seurat_obj, file = "seurat_obj_filtered.Robj")


## EOF
