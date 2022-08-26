# TOOL spatial-transcriptomics-seurat-sctransform.R: "Seurat v4 -Filter spots, normalize with SCTransform and detect high-variance genes in spatial data" (This tool filters out spots with high mitochondrial transcript percentage, indicative of regions where the tissue is damaged. It then normalizes gene expression values using the SCTransform method and detects highly variable genes.)
# INPUT OPTIONAL seurat_spatial_setup.Robj: "Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL seurat_obj_sctransform.Robj
# PARAMETER OPTIONAL mitocutoff: "Filter out spots which have higher mitochondrial transcript percentage" TYPE DECIMAL FROM 0 TO 100 DEFAULT 20 (Filter out spots from regions of damaged tissue. The spots to be kept must have lower percentage of mitochondrial transcripts than this.)
# PARAMETER OPTIONAL hbcutoff: "Filter out spots which have higher hemoglobin transcript percentage" TYPE DECIMAL FROM 0 TO 100 DEFAULT 20 (Filter out spots which have higher percentage of hemoglobin transcripts than this.)
# PARAMETER OPTIONAL num.features: "Number of variable genes to return" TYPE INTEGER DEFAULT 3000 (Number of features to select as top variable features, i.e. how many features returned. For SCTransform, the recommended default is 3000.)
# RUNTIME R-4.1.0-single-cell


# PARAMETER OPTIONAL genes: "Genes to filter" TYPE [NULL: Null, mt: Mitochondrial, Hb: Hemoglobin, both: Both] DEFAULT NULL (You can choose to remove genes based on the top expressed genes plot for example mitochondrial or hemoglobin genes.)

# 2022-07-19 IH

library(Seurat)

# Load the R-Seurat-object (called seurat_obj)
load("seurat_spatial_setup.Robj")

# Subset: remove potential empties, multiplets and broken cells based on parameters
seurat_obj <- PercentageFeatureSet(seurat_obj, "^mt-", col.name = "percent_mito")
seurat_obj <- PercentageFeatureSet(seurat_obj, "^Hb.*-", col.name = "percent_hb")
seurat_obj= seurat_obj[, seurat_obj$percent_mito < mitocutoff & seurat_obj$percent_hb < hbcutoff]

#filter genes
#if (genes == "Hb" | genes == "both") {
#seurat_obj <- seurat_obj[!grepl("^Hb.*-", rownames(seurat_obj)), ]}
#if (genes == "mt" | genes == "both") {
#seurat_obj <- seurat_obj[!grepl("^mt-", rownames(seurat_obj)), ]}

#Sctransform normalizes the data, detects high-variance features, and stores the data in the SCT assay.
seurat_obj <- SCTransform(seurat_obj, assay = "Spatial", variable.features.n = num.features, verbose = FALSE)

# Save the Robj for the next tool
save(seurat_obj, file="seurat_obj_sctransform.Robj")


## EOF
