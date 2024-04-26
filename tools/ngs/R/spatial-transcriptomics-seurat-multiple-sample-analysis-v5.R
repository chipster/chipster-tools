# TOOL spatial-transcriptomics-seurat-multiple-sample-analysis-v5.R: "Seurat v5 -Combine multiple samples, normalization, PCA, clustering and visualization" (This tool performs a combined analysis by merging or integrating multiple data sets across sections. The combined analysis includes similar steps as the single sample analysis which includes normalization, PCA, clustering and visualization. Use this tool instead of the single sample analysis tool when you want to analyze multiple samples at the same time.)
# INPUT samples{...}.Robj: "Samples to combine" TYPE GENERIC
# OUTPUT OPTIONAL seurat_obj_multiple.Robj
# OUTPUT OPTIONAL merged_plot.pdf
# OUTPUT OPTIONAL integrated_plot.pdf
# OUTPUT OPTIONAL PCAloadings.txt
# PARAMETER OPTIONAL num.features: "Number of variable genes to return in SCTransform" TYPE INTEGER DEFAULT 3000 (Number of features to select as top variable features, i.e. how many features returned. For SCTransform, the recommended default is 3000.)
# PARAMETER OPTIONAL method: "Combining method" TYPE [merge: Merge, integration: Integration] DEFAULT merge (User can choose to merge or integrate the samples.)
# PARAMETER OPTIONAL PCstocompute: "Number of PCs to compute" TYPE INTEGER DEFAULT 50 (Number of PCs to compute in PCA.)
# PARAMETER OPTIONAL loadings: "Print loadings in a file" TYPE [TRUE: yes, FALSE: no] DEFAULT FALSE (Print the PC loadings to a txt file.)
# PARAMETER OPTIONAL num.of.genes.loadings: "Number of genes to list in the loadings file" TYPE INTEGER DEFAULT 5 (How many genes to list in the loadings txt file.)
# PARAMETER OPTIONAL dims.reduction: "Dimensions of reduction to use as input" TYPE INTEGER DEFAULT 30 (Number of dimensions of reduction to use for clustering and UMAP. If integration is used, this is also the number of dimensions of reduction to use for integration.)
# PARAMETER OPTIONAL res: "Resolution for granularity for clustering" TYPE DECIMAL DEFAULT 0.8 (Resolution parameter that sets the granularity of the clustering. Increased values lead to greater number of clusters.)
# RUNTIME R-4.2.3-seurat5
# SLOTS 3
# TOOLS_BIN ""

# 2022-07-21 IH
# 2022-10-13 ML nfeatures = 3000 fix
# 2024-03-21 EP Update to Seurat v5 (and add SCTransform to this tool)

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
    # Merging follows spatial analysis vignette (https://satijalab.org/seurat/articles/spatial_vignette.html#10x-visium)
    # where seurat objects are first individually normalized using SCTransform before merging objects 

    # SCTransform
    sctransformed_obj_1 <- SCTransform(seurat_objects[[1]], variable.features.n = num.features, assay = "Spatial", verbose = FALSE)
    sctransformed_obj_2 <- SCTransform(seurat_objects[[2]], variable.features.n = num.features, assay = "Spatial", verbose = FALSE)
    sctransformed_objects <- c(sctransformed_obj_1, sctransformed_obj_2)

    if (length(seurat_objects) >= 3) {
       for (i in 3:length(seurat_objects)) {
            sctransformed_seurat_obj <- SCTransform(seurat_objects[[i]], variable.features.n = num.features, assay = "Spatial", verbose = FALSE)
            sctransformed_objects <- append(sctransformed_objects, sctransformed_seurat_obj)
        }
    }
    
    # Merge multiple slices
    if (length(seurat_objects) == 2) {
        seurat_obj <- merge(sctransformed_objects[[1]], sctransformed_objects[[2]])
    } else {
        seurat_obj <- merge(sctransformed_objects[[1]], y = c(sctransformed_objects[c(2, (length(sctransformed_objects)))]))
    }

    #Set default assay and variable features for merged Seurat object
    DefaultAssay(seurat_obj) <- "SCT"

    variables_list <- c(VariableFeatures(sctransformed_objects[[1]]), VariableFeatures(sctransformed_objects[[2]]))

    if (length(sctransformed_objects) >= 3) {
        for (i in 3:length(sctransformed_objects)) {
            print(variables_list)
            variables_list <- append(variables_list, VariableFeatures(sctransformed_objects[[i]]))
        }
    }
    # Before selecting list of variable features, remove features that are not is scale.data.
    # This is done because without it there will be a 'subsrict out of bounds' error in FindSpatiallyVariableFeatures()
    # as features that are not in both objects are not scaled. This should be taken into account such that the list 
    # of variable features only contain features that are in both objects 
    # More information in issues  https://github.com/satijalab/seurat/issues/3041 and https://github.com/satijalab/seurat/issues/4611
    new_variables_list <- intersect(variables_list,rownames(seurat_obj[['SCT']]$scale.data)) 

    # Select variable features
    VariableFeatures(seurat_obj) <- new_variables_list

    # PCA
    seurat_obj <- RunPCA(seurat_obj, npcs = PCstocompute, assay = "SCT", verbose = FALSE)

    # PCA genes in txt file
    if (loadings == TRUE) {
        sink("PCAloadings.txt")
        print(seurat_obj[["pca"]], dims = 1:PCstocompute, nfeatures = num.of.genes.loadings)
        sink()
    }

    seurat_obj <- FindNeighbors(seurat_obj, dims = 1:dims.reduction, verbose=FALSE)
    seurat_obj <- FindClusters(seurat_obj, resolution = res, verbose = FALSE)
    seurat_obj <- RunUMAP(seurat_obj, dims = 1:dims.reduction,  verbose=FALSE)
    
    # Visualization
    pdf(file = "merged_plot.pdf", , width = 9, height = 12)
    print(DimPlot(seurat_obj, reduction = "umap", group.by = "ident"))
    print(DimPlot(seurat_obj, reduction = "umap", group.by = "orig.ident"))
    print(SpatialDimPlot(seurat_obj, label = TRUE, label.size = 3))

    # Close pdf
    dev.off()
}

if (method == "integration") {
    # Integration follows single-cell analysis vignette (https://satijalab.org/seurat/articles/integration_introduction)
    # where seurat objects are first merged and then individually normalized using SCTransform by using different sample
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
    seurat_obj <- SCTransform(seurat_obj, variable.features.n = num.features, assay = "Spatial", verbose = FALSE)

    # There is an issue with running SCTransform for a merged object as explained by StepahieHowe in issue 
    # https://github.com/satijalab/seurat/issues/8235. This solution does not solve issues with
    # FindSpatiallyVariableFeatures() but it allows to perform integration with scRNA-seq data without
    # any errors. FindSpatiallyVariableFeatures() must still be run for each object individually because
    # otherwise it will result in please provide the same number of observations as spatial locations' error
    # as explained by tingchiafelix
    for (i in length(seurat_objects)) {
        slot(object = seurat_obj@assays$SCT@SCTModel.list[[i]], name="umi.assay")<-"Spatial"
    }

    # You can use this to check that all "umi.assay" slots are called "Spatial"
    print(SCTResults(object=seurat_obj, slot="umi.assay"))

    seurat_obj <- RunPCA(seurat_obj, npcs = PCstocompute, assay = "SCT", verbose = FALSE)

    # PCA genes in txt file
    if (loadings == TRUE) {
        sink("PCAloadings.txt")
        print(seurat_obj[["pca"]], dims = 1:PCstocompute, nfeatures = num.of.genes.loadings)
        sink()
    }

    # Name of new reduction
    new.reduction = "integrated.cca"

    # Integration 
    seurat_obj <- IntegrateLayers(object = seurat_obj, method = "CCAIntegration", normalization.method = "SCT", new.reduction = new.reduction,
    dims = 1:dims.reduction, assay = "SCT", verbose = FALSE)

    seurat_obj <- FindNeighbors(seurat_obj, reduction = new.reduction, dims = 1:dims.reduction, verbose = FALSE)
    seurat_obj <- FindClusters(seurat_obj, resolution = res, verbose = FALSE)
    seurat_obj <- RunUMAP(seurat_obj, reduction = new.reduction, dims = 1:dims.reduction, verbose = FALSE)

    # Visualization
    pdf(file = "integrated_plot.pdf", width = 9, height = 12) 

    print(DimPlot(seurat_obj, reduction = "umap", group.by = "ident"))
    print(DimPlot(seurat_obj, reduction = "umap", group.by = "orig.ident"))
    print(SpatialDimPlot(seurat_obj, label = TRUE, label.size = 3))

    # Close the pdf
    dev.off()
   
}

# Save the combined Robj for the next tool
save(seurat_obj, file = "seurat_obj_multiple.Robj")

# EOF
