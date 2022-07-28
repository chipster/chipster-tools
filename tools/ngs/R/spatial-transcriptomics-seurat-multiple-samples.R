# TOOL spatial-transcriptomics-seurat-multiple-samples.R: "Combine multiple samples" (This tool can be used to combine multiple datasets either by merging or integrating the data across sections. Integration should be done if there are strong batch effects present in the data.)
# INPUT samples{...}.Robj: "Samples to combine" TYPE GENERIC
# OUTPUT OPTIONAL seurat_obj_multiple.Robj
# PARAMETER OPTIONAL method: "Combine the data" TYPE [merge: Merge, integration: Integration] DEFAULT merge (User can choose to merge or integrate the data.)
# RUNTIME R-4.1.0-single-cell
# SLOTS 3

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

# Load the R-Seurat-objects (called seurat_obj -that's why we need to rename them here)
input <- read.table("chipster-inputs.tsv",header = FALSE,sep = "\t")
for (i in 1:nrow(input)) {
    load(input[i,1])
    name <- paste("seurat_obj_", i, sep = "")
    assign(name, seurat_obj)
}
seurat_objects <- as.list(mget(objects(pattern="seurat_obj_"))) 

# Check the number of samples
if (length(seurat_objects) < 2) {
    stop(paste("CHIPSTER-NOTE: ", "It seems you don't have multiple samples. Please check your input files."))
  }

# Merge multiple slices
if (method == "merge") {
    if (length(seurat_objects) == 2) {
        objects_combined <- merge(seurat_objects[[1]], seurat_objects[[2]])
    } else {
        objects_combined <- merge(seurat_objects[[1]], y = c(seurat_objects[c(2,(length(seurat_objects)))]))
    }
    #Set default assay and variable features for merged Seurat object
    DefaultAssay(objects_combined) <- "SCT"
    for (object in seurat_objects) {
        variables_list <- VariableFeatures(object)
    }
    VariableFeatures(objects_combined) <- (c(variables_list))
}

# Integration:
# Code from https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/seurat/seurat_07_spatial.html#Spatial_transcriptomics
if (method == "integration") {
    # need to set maxSize for PrepSCTIntegration to work
    options(future.globals.maxSize = 3000 * 1024^2)  # set allowed size to 3K MiB

    # Select features that are repeatedly variable across datasets for integration
    features <- SelectIntegrationFeatures(seurat_objects, nfeatures = 300, verbose = FALSE)
    # run PrepSCTIntegration to compute the sctransform residuals for all genes in the datasets
    seurat_objects <- PrepSCTIntegration(object.list = seurat_objects, anchor.features = features,
    verbose = FALSE)

    # Identify anchors
    int.anchors <- FindIntegrationAnchors(object.list = seurat_objects, normalization.method = "SCT",
    verbose = FALSE, anchor.features = features)
    # Integrate datasets with the found anchors
    objects_combined <- IntegrateData(anchorset = int.anchors, normalization.method = "SCT",
    verbose = FALSE)

    #remove objects and memory
    rm(int.anchors, seurat_objects)
    gc()
}
#rename combined Robj
seurat_obj <- objects_combined
# Save the combined Robj for the next tool
save(seurat_obj, file="seurat_obj_multiple.Robj") 

#EOF



