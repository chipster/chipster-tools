# TOOL spatial-transcriptomics-seurat-filter-v5.R: "Seurat v5 -Filter spots" (This tool filters out spots based on mitochondrial, ribosomal and hemoglobin transcript percentage.)
# INPUT OPTIONAL seurat_spatial_setup.Robj: "Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL seurat_obj_filtered.Robj
# PARAMETER OPTIONAL mitocutoff: "Filter out spots which have higher mitochondrial transcript percentage" TYPE DECIMAL FROM 0 TO 100 DEFAULT 20 (Filter out spots from regions of damaged tissue. The spots to be kept must have lower percentage of mitochondrial transcripts than this.)
# PARAMETER OPTIONAL minribo: "Filter out cells which have lower ribosomal transcript percentage" TYPE DECIMAL FROM 0 TO 100 DEFAULT 0 (Filter out cells that have lower ribosomal transcript percentage.)
# PARAMETER OPTIONAL hbcutoff: "Filter out spots which have higher hemoglobin transcript percentage" TYPE DECIMAL FROM 0 TO 100 DEFAULT 20 (Filter out spots which have higher percentage of hemoglobin transcripts than this.)
# RUNTIME R-4.2.3-seurat5
# TOOLS_BIN ""


# PARAMETER OPTIONAL genes: "Genes to filter" TYPE [NULL: Null, mt: Mitochondrial, Hb: Hemoglobin, both: Both] DEFAULT NULL (You can choose to remove genes based on the top expressed genes plot for example mitochondrial or hemoglobin genes.)

# 2022-07-19 IH
# 2023-07-12 IH Mitogenes with Mt- and MT- (copied from ML)
# 2023-09-08 IH add ribosomal filtering
# 2024-03-21 EP Update to Seurat v5 (also remove SCTransform and add it to other tool(s))

library(Seurat)
library(Biobase)

print(package.version("Seurat"))
documentVersion("Seurat", package.version("Seurat"))

# Load the R-Seurat-object (called seurat_obj)
load("seurat_spatial_setup.Robj")

print(Version(seurat_obj))

# Subset: remove potential empties, multiplets and broken cells based on parameters
seurat_obj <- subset(seurat_obj, subset = percent.mt <= mitocutoff & percent.hb <= hbcutoff & percent.rb >= minribo)

# Save the Robj for the next tool
save(seurat_obj, file = "seurat_obj_filtered.Robj")


## EOF
