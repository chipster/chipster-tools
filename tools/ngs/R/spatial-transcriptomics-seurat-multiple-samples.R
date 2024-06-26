# TOOL spatial-transcriptomics-seurat-multiple-samples.R: "Seurat v4 -Combine multiple samples" (This tool can be used to combine multiple datasets either by merging or integrating the data across sections. Integration should be done if there are strong batch effects present in the data.)
# INPUT samples{...}.Robj: "Samples to combine" TYPE GENERIC
# OUTPUT OPTIONAL seurat_obj_multiple.Robj
# OUTPUT OPTIONAL logfile.txt
# PARAMETER OPTIONAL method: "Combining method" TYPE [merge: Merge, integration: Integration] DEFAULT merge (User can choose to merge or integrate the samples.)
# PARAMETER OPTIONAL number_of_var_feats: "Number of variable genes in combined object " TYPE INTEGER DEFAULT 3000 (Number of variable genes in merged/integrated object. If integration is used, this is also the number of genes used for integration. This value should be less than or equal to the number of variable genes returned in SCTransform. This tool selects a list of highest variable genes across all data sets based on the user selected number. PCA is only run on these variable genes. In case this number of variable genes is set too high, there may not be as many variable genes across all data sets. The tool will still run, but a log file will be produced to inform the user.)
# RUNTIME R-4.2.3-single-cell
# SLOTS 3
# TOOLS_BIN ""

# 2022-07-21 IH
# 2022-10-13 ML nfeatures = 3000 fix

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

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

# Merge multiple slices
if (method == "merge") {
    if (length(seurat_objects) == 2) {
        objects_combined <- merge(seurat_objects[[1]], seurat_objects[[2]])
    } else {
        objects_combined <- merge(seurat_objects[[1]], y = c(seurat_objects[c(2, (length(seurat_objects)))]))
    }
    # Set default assay and variable features for merged Seurat object
    DefaultAssay(seurat_obj) <- "SCT"

    # Select variable features using SelectIntegrationFeatures() as explained in https://github.com/satijalab/seurat/issues/5761
    for (seurat_obj in seurat_objects) {
        if (number_of_var_feats > length(VariableFeatures(seurat_obj))) {
            stop(paste("CHIPSTER-NOTE: ", "The selected number of variable genes in the combined object should be smaller than or equal to the selected number of returned variable genes in SCTransform in the previous tool."))
        }
    }
    variable_feats <- SelectIntegrationFeatures(seurat_objects, nfeatures = number_of_var_feats, verbose = FALSE)

    # Variable genes that are not in all data sets (scale.data) need to be removed from the list before downstream analyses.
    # This is done because PCA and FindSpatiallyVariableFeatures() exclude genes not in scale.data from analysis.
    feats_in_all_objects <- intersect(rownames(seurat_objects[[1]]$SCT@scale.data), rownames(seurat_objects[[2]]$SCT@scale.data))  
    variable_feats_in_all_objects <- intersect(variable_feats, feats_in_all_objects)
    if (length(variable_feats_in_all_objects) < length(variable_feats)) {
        sink("logfile.txt")
        warning_message <- paste("
        You have selected more variable genes for the combined object than have been deemed variable across
        all data sets. PCA in the next tool is only run with", length(variable_feats_in_all_objects), "variable genes. If you want to rerun the tool, 
        you can try to select a higher number of returned variable genes in SCTransform in the previous tool to see 
        if the number of deemed variable genes across all data sets increases. You can also run the tool again 
        with less than", length(variable_feats_in_all_objects), "selected variable genes for the combined object.")
        cat(warning_message)

        # Close the log file
        sink()
    }
    VariableFeatures(objects_combined) <- variable_feats_in_all_objects
}

# Integration:
# Code mainly from https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/seurat/seurat_07_spatial.html#Spatial_transcriptomics
if (method == "integration") {
    # need to set maxSize for PrepSCTIntegration to work
    options(future.globals.maxSize = 2000 * 1024^2) # set allowed size to 3K MiB

    # Select variable features using SelectIntegrationFeatures() as explained in https://github.com/satijalab/seurat/issues/5761
    for (seurat_obj in seurat_objects) {
        if (number_of_var_feats > length(VariableFeatures(seurat_obj))) {
            stop(paste("CHIPSTER-NOTE: ", "The selected number of variable genes in the combined object should be smaller than or equal to the selected number of returned variable genes in SCTransform in the previous tool."))
        }
    }
    variable_feats <- SelectIntegrationFeatures(seurat_objects, nfeatures = number_of_var_feats, verbose = FALSE)

    # Select features that are repeatedly variable across datasets for integration
    # Variable genes that are not in all data sets (scale.data) need to be removed from the list before downstream analyses.
    # This is done because PCA and FindSpatiallyVariableFeatures() exclude genes not in scale.data from analysis. 
    feats_in_all_objects <- intersect(rownames(seurat_objects[[1]]$SCT@scale.data), rownames(seurat_objects[[2]]$SCT@scale.data))
    variable_feats_in_all_objects <- intersect(variable_feats, feats_in_all_objects)
    if (length(variable_feats_in_all_objects) < length(variable_feats)) {
        sink("logfile.txt")
        warning_message <- paste("
        You have selected more variable genes for the combined object than have been deemed variable across
        all data sets. PCA in the next tool is only run with", length(variable_feats_in_all_objects), "variable genes. If you want to rerun the tool, 
        you can try to select a higher number of returned variable genes in SCTransform in the previous tool to see 
        if the number of deemed variable genes across all data sets increases. You can also run the tool again 
        with less than", length(variable_feats_in_all_objects), "selected variable genes for the combined object.")
        cat(warning_message)

        # Close the log file
        sink()
    }

    # Run PrepSCTIntegration to compute the sctransform residuals for all genes in the datasets
    seurat_objects <- PrepSCTIntegration(
        object.list = seurat_objects, anchor.features = variable_feats_in_all_objects,
        verbose = FALSE
    )
    # Identify anchors
    int.anchors <- FindIntegrationAnchors(
        object.list = seurat_objects, normalization.method = "SCT",
        verbose = FALSE, anchor.features = variable_feats_in_all_objects
    )
    # Integrate datasets with the found anchors
    objects_combined <- IntegrateData(
        anchorset = int.anchors, normalization.method = "SCT",
        verbose = FALSE
    )
    # Remove objects and memory
    rm(int.anchors, seurat_objects)
    gc()
}
# Rename combined Robj
seurat_obj <- objects_combined

# Save the combined Robj for the next tool
save(seurat_obj, file = "seurat_obj_multiple.Robj")

# EOF
