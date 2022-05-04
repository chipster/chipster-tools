# TOOL single-cell-seurat-combine-multiple-samples.R: "Seurat v4 -Combine multiple samples" (This tool can be used to integrate data and combine multiple Seurat objects for later joined analysis. The samples \/R-objects need to be named when created in the Seurat Setup tool.) 
# INPUT samples{...}.Robj: "Samples to combine and align" TYPE GENERIC
# OUTPUT OPTIONAL CCAplot.pdf
# OUTPUT seurat_obj_combined.Robj
# PARAMETER OPTIONAL normalisation.method: "Normalisation method used previously" TYPE [normal:"Global scaling normalization", sctransform:SCTransform] DEFAULT normal (Which normalisation method was used in preprocessing, Global scaling normalization \(default, NormalizeData function used\) or SCTransform.)
# PARAMETER OPTIONAL anchor.identification.method: "Anchor identification method" TYPE [normal:"CCA", rpcamethod:RPCA] DEFAULT normal (Which anchor identification method to use. By default, canonical correlation analysis CCA is used, but user can also decide to use the faster and more conservative reciprocal PCA approach. Check from the manual in which cases this option is recommended.)
# PARAMETER OPTIONAL CCstocompute: "Number of CCs to use in the neighbor search" TYPE INTEGER DEFAULT 20 (Which dimensions to use from the CCA to specify the neighbor search space. The neighbors are used to determine the anchors for the alignment.)
# PARAMETER OPTIONAL PCstocompute: "Number of PCs to use in the anchor weighting" TYPE INTEGER DEFAULT 20 (Number of PCs to use in the anchor weighting procedure. The anchors and their weights are used to compute the correction vectors, which allow the datasets to be integrated.)
# RUNTIME R-4.1.0-single-cell
# SLOTS 3


# 2021-12-30 ML
# 2022-02-17 EK increased slots to 4
# 2022-04-19 ML increased slots to 5
# 2022-05-04 ML add RPCA option for anchor identification



library(Seurat)

# Load the R-Seurat-objects (called seurat_obj -that's why we need to rename them here)
input.names <- read.table("chipster-inputs.tsv",header = FALSE,sep = "\t")
for (i in 1:nrow(input.names)) {
    # unzipIfGZipFile(input.names[i,1])
    load(input.names[i,1])
    name.of.obj <- paste("seurat_obj_", i, sep = "")
    assign(name.of.obj, seurat_obj)
}

# seurat.objects.list <- objects(pattern="seurat_obj")
seurat.objects.list <- as.list(mget(objects(pattern="seurat_obj"))) 

# For testing, please ignore:
# OUTPUT OPTIONAL seurat_obj_list.Robj
# save(seurat.objects.list, file="seurat_obj_list.Robj")

# Select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = seurat.objects.list)


# Perform integration: 

# When data is normalised with NormalizeData:
if (normalisation.method == "normal"){
    # 1. identify anchors using the FindIntegrationAnchors function
    # CCA:
    if (anchor.identification.method == "normal"){
    data.anchors <- FindIntegrationAnchors(object.list = seurat.objects.list, dims = 1:CCstocompute, anchor.features = features) # dims = Which dimensions to use from the CCA to specify the neighbor search space
    # RPCA:
    } else{
        # When using RPCA, run PCA on each dataset using these features
        seurat.objects.list <- lapply(X = seurat.objects.list, FUN = function(x) {
            x <- ScaleData(x, features = features, verbose = FALSE)
            x <- RunPCA(x, features = features, verbose = FALSE)
        })
        
        data.anchors <- FindIntegrationAnchors(object.list = seurat.objects.list, dims = 1:CCstocompute, anchor.features = features, reduction = "rpca") # dims = Which dimensions to use from the CCA to specify the neighbor search space
    }
    # 2. use these anchors to integrate the two datasets together with IntegrateData.
    data.combined <- IntegrateData(anchorset = data.anchors, dims = 1:PCstocompute) # dims = Number of PCs to use in the weighting procedure

    DefaultAssay(data.combined) <- "integrated"

    # Note: these steps are now done twice?
    data.combined <- ScaleData(data.combined, verbose = FALSE)  
    data.combined <- RunPCA(data.combined, npcs = 30, verbose = FALSE)

# When data is normalised with SCTransform:
} else{
    # When SCTransform was used to normalise the data, do a prep step:
    seurat.objects.list <- PrepSCTIntegration(object.list = seurat.objects.list, anchor.features = features)
    # 1. identify anchors using the FindIntegrationAnchors function
     # CCA:
    if (anchor.identification.method == "normal"){
    data.anchors <- FindIntegrationAnchors(object.list = seurat.objects.list, dims = 1:CCstocompute, anchor.features = features, normalization.method = "SCT") # dims = Which dimensions to use from the CCA to specify the neighbor search space
    # RPCA:
    } else{
        # When using RPCA, run PCA on each dataset using these features
        seurat.objects.list <- lapply(X = seurat.objects.list, FUN = function(x) {
            x <- ScaleData(x, features = features, verbose = FALSE)
            x <- RunPCA(x, features = features, verbose = FALSE)
        })
        
        data.anchors <- FindIntegrationAnchors(object.list = seurat.objects.list, dims = 1:CCstocompute, anchor.features = features, normalization.method = "SCT", reduction = "rpca") # dims = Which dimensions to use from the CCA to specify the neighbor search space
    }
    
    
    # 2. use these anchors to integrate the two datasets together with IntegrateData.
    data.combined <- IntegrateData(anchorset = data.anchors, dims = 1:PCstocompute, normalization.method = "SCT") # dims = Number of PCs to use in the weighting procedure

    DefaultAssay(data.combined) <- "integrated"

    # Note: Skip ScaleData when using SCTransform
    data.combined <- ScaleData(data.combined, verbose = FALSE)  
    data.combined <- RunPCA(data.combined, npcs = 30, verbose = FALSE)
}


# Save the Robj for the next tool
save(data.combined, file="seurat_obj_combined.Robj")

## EOF



