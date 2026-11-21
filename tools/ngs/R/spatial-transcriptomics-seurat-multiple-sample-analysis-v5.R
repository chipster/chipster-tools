# TOOL spatial-transcriptomics-seurat-multiple-sample-analysis-v5.R: "Seurat v5 -Combine multiple samples, normalization and PCA" (This tool combines multiple Seurat objects into one object, normalizes each sample with SCTransform, selects the variable genes across the samples and runs PCA. Normalized expression values are the same for both combining methods. If integration is selected, the samples are also integrated with CCA, which adds an integrated embedding \(integrated.cca\) that later tools use instead of PCA. Integration does not change the expression values.)
# INPUT samples{...}.Robj: "Samples to combine" TYPE GENERIC
# OUTPUT OPTIONAL seurat_obj_multiple.Robj
# OUTPUT OPTIONAL PCAloadings.txt
# OUTPUT OPTIONAL PCAplots.pdf
# PARAMETER OPTIONAL method: "Combining method" TYPE [merge: Merge, integration: Integration] DEFAULT merge (Choose whether to only merge the samples or to also integrate them. Integrate replicates or conditions of the same tissue. Don't integrate sections of different anatomical regions, because integration can then remove real biological differences.)
# PARAMETER OPTIONAL number_of_var_feats: "Number of variable genes" TYPE INTEGER DEFAULT 3000 (Number of variable genes selected across the samples. PCA is run on these genes, and if integration is used, they are also used for integration.)
# PARAMETER OPTIONAL PCstocompute: "Number of PCs to compute" TYPE INTEGER DEFAULT 50 (Number of PCs to compute in PCA. If integration is used, this is also the number of PCs used for integration.)
# PARAMETER OPTIONAL loadings: "Print loadings in a file" TYPE [TRUE: yes, FALSE: no] DEFAULT FALSE (Print the PC loadings in a txt file.)
# PARAMETER OPTIONAL num.of.genes.loadings: "Number of genes to list in the loadings file" TYPE INTEGER DEFAULT 5 (Number of genes to list in the loadings txt file.)
# RUNTIME R-4.2.3-seurat5
# SLOTS 3
# TOOLS_BIN ""

# 2022-07-21 IH
# 2022-10-13 ML nfeatures = 3000 fix
# 2024-03-21 EP Update to Seurat v5 (and add SCTransform to this tool)
# 2025-04-16 ML Fix >3 samples bugs
# 2026-11-09 ML Same normalization for merge and integration: one SCTransform on the merged, layered object

library(Seurat)
library(gplots)
library(ggplot2)
library(patchwork)
library(dplyr)
library(Biobase)
require(cowplot)
library(Matrix)

print(package.version("Seurat"))
documentVersion("Seurat", package.version("Seurat"))

# Load the R-Seurat-objects (called seurat_obj -that's why we need to rename them here)
input <- read.table("chipster-inputs.tsv", header = FALSE, sep = "\t")
for (i in 1:nrow(input)) {
    load(input[i, 1])
    name <- paste("seurat_obj_", i, sep = "")
    assign(name, seurat_obj)
}
seurat_objects <- as.list(mget(objects(pattern = "seurat_obj_")))

# Check the number of samples
if (length(seurat_objects) < 2) {
    stop(paste("CHIPSTER-NOTE: ", "It seems you don't have multiple samples. Please check your input files."))
}

# Merge the unnormalized samples. In Seurat v5 the merged Spatial assay keeps each sample in its own layer
# (counts.1, counts.2, ...), so the SCTransform below still fits a separate model for each sample, like in
# https://satijalab.org/seurat/articles/integration_introduction. No need to split the object first.
seurat_obj <- merge(seurat_objects[[1]], y = seurat_objects[-1])

# SCTransform. Print each different warning once to the job output without stopping SCTransform or running it again.
printed_warnings <- character(0)
seurat_obj <- withCallingHandlers(
    SCTransform(seurat_obj, variable.features.n = number_of_var_feats, assay = "Spatial", verbose = FALSE),
    warning = function(w) {
        if (!conditionMessage(w) %in% printed_warnings) {
            print(w)
            printed_warnings <<- c(printed_warnings, conditionMessage(w))
        }
        invokeRestart("muffleWarning")
    }
)

# There is an issue with running SCTransform for a merged object as explained by StephanieHowe in issue
# https://github.com/satijalab/seurat/issues/8235: the models point to the layer names instead of the
# Spatial assay. Fixing it allows integration with scRNA-seq data and PrepSCTFindMarkers() without errors.
# FindSpatiallyVariableFeatures() must still be run for each image separately, otherwise it fails with
# 'please provide the same number of observations as spatial locations'.
for (i in seq_along(seurat_obj[["SCT"]]@SCTModel.list)) {
    slot(object = seurat_obj[["SCT"]]@SCTModel.list[[i]], name = "umi.assay") <- "Spatial"
}
print(SCTResults(object = seurat_obj, slot = "umi.assay"))

DefaultAssay(seurat_obj) <- "SCT"

# SCTransform selects the variable genes across all samples. It computes the residuals (scale.data) for all of
# them, so all variable genes can be used in PCA.
variable_feats <- head(VariableFeatures(seurat_obj), number_of_var_feats)
variable_feats_scaled <- intersect(variable_feats, rownames(seurat_obj[["SCT"]]$scale.data))
if (length(variable_feats_scaled) < length(variable_feats)) {
    print(paste(
        "Only", length(variable_feats_scaled), "of the", length(variable_feats),
        "variable genes have scaled values, PCA is run with these."
    ))
}
VariableFeatures(seurat_obj) <- variable_feats_scaled

# PCA
seurat_obj <- RunPCA(seurat_obj, npcs = PCstocompute, assay = "SCT", verbose = FALSE)

# PCA genes in txt file
if (loadings == TRUE) {
    sink("PCAloadings.txt")
    print(seurat_obj[["pca"]], dims = 1:PCstocompute, nfeatures = num.of.genes.loadings)
    sink()
}

# PDF plots
pdf(file = "PCAplots.pdf", width = 9, height = 12)
print(ElbowPlot(seurat_obj, ndims = PCstocompute) + ggtitle("Amount of variation in the data explained by each PC"))
print(VizDimLoadings(seurat_obj, dims = 1:1, reduction = "pca") + ggtitle("Top 30 genes associated with PC 1"))
print(VizDimLoadings(seurat_obj, dims = 2:2, reduction = "pca") + ggtitle("Top 30 genes associated with PC 2"))
dev.off() # close the pdf

# INTEGRATION
# In Seurat v5, IntegrateLayers() only adds the integrated embedding (integrated.cca). It doesn't change the
# counts, data or scale.data of the SCT assay. The clustering tool uses integrated.cca when the object has it.
if (method == "integration") {
    seurat_obj <- IntegrateLayers(
        object = seurat_obj,
        method = CCAIntegration,
        normalization.method = "SCT",
        new.reduction = "integrated.cca",
        dims = 1:PCstocompute,
        assay = "SCT",
        features = variable_feats_scaled,
        verbose = FALSE
    )
}

# Save the combined Robj for the next tool
save(seurat_obj, file = "seurat_obj_multiple.Robj")

# EOF
