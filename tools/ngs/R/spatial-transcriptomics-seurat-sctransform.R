# TOOL spatial-transcriptomics-seurat-sctransform.R: "Seurat v4 -SCTransform: Filter cells, normalize, regress and detect high-variance features in spatial data" (This tool filters out dead cells, empties and doublets. It then normalizes gene expression values using the SCTransform method, detects highly variable genes, scales the data and regresses out unwanted variation based on the number of UMIs and mitochondrial transcript percentage. You can also choose to regress out variation due to cell cycle heterogeneity.)
# INPUT OPTIONAL seurat_spatial_setup.Robj: "Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL seurat_obj_sctransform.Robj
# PARAMETER OPTIONAL mitocutoff: "Filter out cells which have higher mitochondrial transcript percentage" TYPE DECIMAL FROM 0 TO 100 DEFAULT 20 (Filter out dead cells. The cells to be kept must have lower percentage of mitochondrial transcripts than this.)
# PARAMETER OPTIONAL hbcutoff: "Filter out cells which have higher hemoglobin transcript percentage" TYPE DECIMAL FROM 0 TO 100 DEFAULT 20 (Filter out dead cells. The cells to be kept must have lower percentage of hemoglobin transcripts than this.)
# PARAMETER OPTIONAL num.features: "Number of variable genes to return" TYPE INTEGER DEFAULT 3000 (Number of features to select as top variable features, i.e. how many features returned. For SCTransform, the recommended default is 3000.)
# RUNTIME R-4.1.0-single-cell


# PARAMETER OPTIONAL genes: "Genes to filter" TYPE [NULL: Null, mt: Mitochondrial, Hb: Hemoglobin, both: Both] DEFAULT NULL (You can choose to remove genes based on the top expressed genes plot for example mitochondrial or hemoglobin genes.)


library(Seurat)
library(SeuratData)

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
