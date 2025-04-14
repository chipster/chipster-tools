# TOOL spatial-transcriptomics-seurat-multiple-sample-analysis-v5.R: "Seurat v5 -Combine multiple samples, normalization and PCA" (This tool combines multiple Seurat objects and performs normalization and PCA for multiple samples at the same time. This tool first normalizes each sample with SCTransform and detects highly variable genes. Then, the normalized samples are combined by merging the samples. After combining the samples, variable genes for the combined Seurat object are selected. Then, PCA is performed on the merged object with the selected variable genes. Finally, if integration is selected, samples are integrated using CCA as the anchor identification method.)
# INPUT samples{...}.Robj: "Samples to combine" TYPE GENERIC
# OUTPUT OPTIONAL seurat_obj_multiple.Robj
# OUTPUT OPTIONAL PCAloadings.txt
# OUTPUT OPTIONAL logfile.txt
# OUTPUT OPTIONAL PCAplots.pdf
# PARAMETER OPTIONAL num.features: "Number of variable genes to return in SCTransform" TYPE INTEGER DEFAULT 3000 (Number of variable genes to return in SCTransform.)
# PARAMETER OPTIONAL method: "Combining method" TYPE [merge: Merge, integration: Integration] DEFAULT merge (Choose whether to merge or integrate the samples.)
# PARAMETER OPTIONAL number_of_var_feats: "Number of variable genes in combined object " TYPE INTEGER DEFAULT 3000 (Number of variable genes in merged/integrated object. If integration is used, this is also the number of genes used for integration. This value should be less than or equal to the number of variable genes returned in SCTransform. This tool selects a list of highest variable genes across all data sets based on the user selected number. PCA is only run on these variable genes. In case this number of variable genes is set too high, there may not be as many variable genes across all data sets after SCTransform. The tool will still run, but a log file will be produced to inform the user.)
# PARAMETER OPTIONAL PCstocompute: "Number of PCs to compute" TYPE INTEGER DEFAULT 50 (Number of PCs to compute in PCA. If integration is used, this is also the number of PCs used for integration.)
# PARAMETER OPTIONAL loadings: "Print loadings in a file" TYPE [TRUE: yes, FALSE: no] DEFAULT FALSE (Print the PC loadings in a txt file.)
# PARAMETER OPTIONAL num.of.genes.loadings: "Number of genes to list in the loadings file" TYPE INTEGER DEFAULT 5 (Number of genes to list in the loadings txt file.)
# RUNTIME R-4.2.3-seurat5
# SLOTS 3
# TOOLS_BIN ""

# 2022-07-21 IH
# 2022-10-13 ML nfeatures = 3000 fix
# 2024-03-21 EP Update to Seurat v5 (and add SCTransform to this tool)
# 2025-04-14 ML Bug fix in >= 3 sample case

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

if (method == "merge") {
    # Merging mainly follows spatial analysis vignette (https://satijalab.org/seurat/articles/spatial_vignette.html#10x-visium)
    # where seurat objects are first individually normalized using SCTransform before merging objects

    # Run SCTransform and print warnings after SCTransform() but suppress them from Chipster note
    sctransformed_obj_1 <- tryCatch(
        {
            SCTransform(seurat_objects[[1]], variable.features.n = num.features, assay = "Spatial", verbose = FALSE)
        },
        warning = function(w) {
            print(w)
            return(suppressWarnings(SCTransform(seurat_objects[[1]], variable.features.n = num.features, assay = "Spatial", verbose = FALSE)))
        }
    )
    sctransformed_obj_2 <- tryCatch(
        {
            SCTransform(seurat_objects[[2]], variable.features.n = num.features, assay = "Spatial", verbose = FALSE)
        },
        warning = function(w) {
            print(w)
            return(suppressWarnings(SCTransform(seurat_objects[[2]], variable.features.n = num.features, assay = "Spatial", verbose = FALSE)))
        }
    )
    sctransformed_objects <- c(sctransformed_obj_1, sctransformed_obj_2)

    if (length(seurat_objects) >= 3) {
        for (i in 3:length(seurat_objects)) {
            sctransformed_obj <- tryCatch(
                {
                    SCTransform(seurat_objects[[i]], variable.features.n = num.features, assay = "Spatial", verbose = FALSE)
                },
                warning = function(w) {
                    print(w)
                    return(suppressWarnings(SCTransform(seurat_objects[[i]], variable.features.n = num.features, assay = "Spatial", verbose = FALSE)))
                }
            )
            sctransformed_objects <- append(sctransformed_objects, sctransformed_obj)
        }
    }

    # Merge multiple slices
    if (length(seurat_objects) == 2) {
        seurat_obj <- merge(sctransformed_objects[[1]], sctransformed_objects[[2]])
    } else {
        seurat_obj <- merge(sctransformed_objects[[1]], y = c(sctransformed_objects[c(2, (length(sctransformed_objects)))]))
    }

    # Set default assay and variable features for merged Seurat object
    DefaultAssay(seurat_obj) <- "SCT"

    # Select variable features using SelectIntegrationFeatures() as explained in https://github.com/satijalab/seurat/issues/5761
    if (number_of_var_feats > num.features) {
        stop(paste("CHIPSTER-NOTE: ", "The selected number of variable genes in the combined object should be smaller than or equal to the selected number of returned variable genes in SCTransform."))
    } else {
        variable_feats <- SelectIntegrationFeatures(sctransformed_objects, nfeatures = number_of_var_feats, verbose = FALSE)
    }

    # Variable genes that are not in all data sets (scale.data) need to be removed from the list before downstream analyses.
    # This is done because without it there will be a 'subsrict out of bounds' error in FindSpatiallyVariableFeatures()
    # as genes that are not in all objects are not scaled. In addition, PCA excludes genes not in scale.data from analysis.
    # More information in issues  https://github.com/satijalab/seurat/issues/3041 and https://github.com/satijalab/seurat/issues/4611
    variable_feats_in_all_objects <- intersect(variable_feats, rownames(seurat_obj[["SCT"]]$scale.data))
    if (length(variable_feats_in_all_objects) < length(variable_feats)) {
        sink("logfile.txt")
        warning_message <- paste("
        You have selected more variable genes for the combined object than have been deemed variable across
        all data sets after SCTransform. PCA is only run with", length(variable_feats_in_all_objects), "variable genes. If you want to rerun the tool,
        you can try to select a higher number of returned variable genes in SCTransform to see if the number
        of deemed variable genes across all data sets increases. You can also run the tool again with less
        than", length(variable_feats_in_all_objects), "selected variable genes for the combined object.")
        cat(warning_message)

        # Close the log file
        sink()
    }
    VariableFeatures(seurat_obj) <- variable_feats_in_all_objects

    # PCA
    seurat_obj <- RunPCA(seurat_obj, npcs = PCstocompute, assay = "SCT", verbose = FALSE)

    # PCA genes in txt file
    if (loadings == TRUE) {
        sink("PCAloadings.txt")
        print(seurat_obj[["pca"]], dims = 1:PCstocompute, nfeatures = num.of.genes.loadings)
        sink()
    }
    # PDF plots
    pdf(file = "PCAplots.pdf", , width = 9, height = 12)
    print(ElbowPlot(seurat_obj, ndims = PCstocompute) + ggtitle("Amount of variation in the data explained by each PC"))
    print(VizDimLoadings(seurat_obj, dims = 1:1, reduction = "pca") + ggtitle("Top 30 genes associated with PC 1"))
    print(VizDimLoadings(seurat_obj, dims = 2:2, reduction = "pca") + ggtitle("Top 30 genes associated with PC 2"))
    dev.off() # close the pdf
}

if (method == "integration") {
    # Integration follows single-cell analysis vignette (https://satijalab.org/seurat/articles/integration_introduction)
    # where a combined object is used for SCTransform, but samples are individually normalized using SCTransform by using different sample
    # layers stored in the object. Seurat object needs to be split by the samples before SCTransform so that SCTransform
    # is run for each sample individually. However, splitting the object as shown in the vignette is not needed after merging
    # multiple unnormalized Seurat objects because objects are already split by their sample (=orig.ident) after merge

    # Merge multiple slices
    if (length(seurat_objects) == 2) {
        seurat_obj <- merge(seurat_objects[[1]], seurat_objects[[2]])
    } else {
        seurat_obj <- merge(seurat_objects[[1]], y = c(seurat_objects[c(2, (length(seurat_objects)))]))
    }

    # SCTransform (no need to split object before SCTransform)
    # This will also set variable genes in all objects (different compared to merge)
    # Print warnings after SCTransform() but suppress them from Chipster note
    seurat_obj <- tryCatch(
        {
            SCTransform(seurat_obj, variable.features.n = num.features, assay = "Spatial", verbose = FALSE)
        },
        warning = function(w) {
            print(w)
            return(suppressWarnings(SCTransform(seurat_obj, variable.features.n = num.features, assay = "Spatial", verbose = FALSE)))
        }
    )

    # There is an issue with running SCTransform for a merged object as explained by StepahieHowe in issue
    # https://github.com/satijalab/seurat/issues/8235. This solution does not solve issues with
    # FindSpatiallyVariableFeatures() but it allows to perform integration with scRNA-seq data without
    # any errors. FindSpatiallyVariableFeatures() must still be run for each object individually because
    # otherwise it will result in 'please provide the same number of observations as spatial locations' error
    # as explained by tingchiafelix
    for (i in length(seurat_objects)) {
        slot(object = seurat_obj@assays$SCT@SCTModel.list[[i]], name = "umi.assay") <- "Spatial"
    }

    # You can use this to check that all "umi.assay" slots are called "Spatial"
    print(SCTResults(object = seurat_obj, slot = "umi.assay"))

    # Set default assay and variable features before integration
    DefaultAssay(seurat_obj) <- "SCT"

    # Select variable features for combined object from combined object after SCTransform
    if (number_of_var_feats > num.features) {
        stop(paste("CHIPSTER-NOTE: ", "The selected number of variable genes in the combined object should be smaller than or equal to the selected number of returned variable genes in SCTransform."))
    } else {
        variable_feats <- VariableFeatures(seurat_obj)[1:number_of_var_feats]
    }

    # Variable genes that are not in all data sets (scale.data) need to be removed from the list before downstream analyses.
    # This is done because without it there will be a 'subsrict out of bounds' error in FindSpatiallyVariableFeatures()
    # as genes that are not in all objects are not scaled. In addition, PCA excludes genes not in scale.data from analysis.
    # More information in issues  https://github.com/satijalab/seurat/issues/3041 and https://github.com/satijalab/seurat/issues/4611
    variable_feats_in_all_objects <- intersect(variable_feats, rownames(seurat_obj[["SCT"]]$scale.data))
    if (length(variable_feats_in_all_objects) < length(variable_feats)) {
        sink("logfile.txt")
        warning_message <- paste("
        You have selected more variable genes for the combined object than have been deemed variable across
        all data sets after SCTransform. PCA is only run with", length(variable_feats_in_all_objects), "variable genes. If you want to rerun the tool,
        you can try to select a higher number of returned variable genes in SCTransform to see if the number
        of deemed variable genes across all data sets increases. You can also run the tool again with less
        than", length(variable_feats_in_all_objects), "selected variable genes for the combined object.")
        cat(warning_message)
        # Close the log file
        sink()
    }
    VariableFeatures(seurat_obj) <- variable_feats_in_all_objects

    seurat_obj <- RunPCA(seurat_obj, npcs = PCstocompute, assay = "SCT", verbose = FALSE)

    # PCA genes in txt file
    if (loadings == TRUE) {
        sink("PCAloadings.txt")
        print(seurat_obj[["pca"]], dims = 1:PCstocompute, nfeatures = num.of.genes.loadings)
        sink()
    }

    # PDF plots
    pdf(file = "PCAplots.pdf", , width = 9, height = 12)
    print(ElbowPlot(seurat_obj, ndims = PCstocompute) + ggtitle("Amount of variation in the data explained by each PC"))
    print(VizDimLoadings(seurat_obj, dims = 1:1, reduction = "pca") + ggtitle("Top 30 genes associated with PC 1"))
    print(VizDimLoadings(seurat_obj, dims = 2:2, reduction = "pca") + ggtitle("Top 30 genes associated with PC 2"))
    dev.off() # close the pdf

    # Name of new reduction
    new.reduction <- "integrated.cca"

    # Integration
    seurat_obj <- IntegrateLayers(
        object = seurat_obj, method = "CCAIntegration", normalization.method = "SCT", new.reduction = new.reduction,
        dims = 1:PCstocompute, assay = "SCT", features = variable_feats_in_all_objects, verbose = FALSE
    )
}

# Save the combined Robj for the next tool
save(seurat_obj, file = "seurat_obj_multiple.Robj")

# EOF
