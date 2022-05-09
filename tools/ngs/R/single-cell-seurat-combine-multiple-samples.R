# TOOL single-cell-seurat-combine-multiple-samples.R: "Seurat v4 -Combine multiple samples" (This tool can be used to integrate data and combine multiple Seurat objects for later joined analysis. The samples \/R-objects need to be named when created in the Seurat Setup tool.) 
# INPUT samples{...}.Robj: "Samples to combine and align" TYPE GENERIC
# OUTPUT OPTIONAL CCAplot.pdf
# OUTPUT seurat_obj_combined.Robj
# PARAMETER OPTIONAL normalisation.method: "Normalisation method used previously" TYPE [LogNormalize:"Global scaling normalization", SCT:"SCTransform"] DEFAULT LogNormalize (Which normalisation method was used in preprocessing, Global scaling normalization \(default, NormalizeData function used\) or SCTransform.)
# PARAMETER OPTIONAL anchor.identification.method: "Anchor identification method" TYPE [cca:CCA, rpca:RPCA] DEFAULT cca (Which anchor identification method to use. By default, canonical correlation analysis CCA is used, but user can also decide to use the faster and more conservative reciprocal PCA approach. Check from the manual in which cases this option is recommended.)
# PARAMETER OPTIONAL CCstocompute: "Number of CCs to use in the neighbor search" TYPE INTEGER DEFAULT 20 (Which dimensions to use from the CCA to specify the neighbor search space. The neighbors are used to determine the anchors for the alignment.)
# PARAMETER OPTIONAL PCstocompute: "Number of PCs to use in the anchor weighting" TYPE INTEGER DEFAULT 20 (Number of PCs to use in the anchor weighting procedure. The anchors and their weights are used to compute the correction vectors, which allow the datasets to be integrated.)
# PARAMETER OPTIONAL ref.sample.names: "Samples to use as references" TYPE STRING DEFAULT "No references selected" (Names of the sample or samples you wish to use as references in integration, separated by comma. If you are integrating several large datasets, the tool might run out of memory. Choosing to use only some of them as references makes the integration more memory efficient and faster. Please note that the sample names here are case sensitive, so check how you typed the names of the samples when running the setup tool.)
# RUNTIME R-4.1.0-single-cell
# SLOTS 3


# 2021-12-30 ML
# 2022-02-17 EK increased slots to 4
# 2022-04-19 ML increased slots to 5
# 2022-05-04 ML add RPCA option for anchor identification
# 2022-05-05 ML Rewrite the code, add option to use only part of samples as references



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

# When SCTransform was used to normalise the data, do a prep step:
if (normalisation.method == "SCT"){
    seurat.objects.list <- PrepSCTIntegration(object.list = seurat.objects.list, anchor.features = features)
}

# When using RPCA, need to run PCA on each dataset using these features:
if (anchor.identification.method == "rpca"){
        seurat.objects.list <- lapply(X = seurat.objects.list, FUN = function(x) {
            x <- ScaleData(x, features = features, verbose = FALSE)
            x <- RunPCA(x, features = features, verbose = FALSE)
        })
}

# When using only the user listed samples as references:
if (ref.sample.names != "No references selected"){
    ref.samples.name.list <- unlist(strsplit(ref.sample.names, ", "))
    # Go through the samples = R-objects in the list to see which ones are the reference samples.
    ref.sample.numbers <- vector()
    for (i in 1:length(seurat.objects.list)) {
     # Check, if the (first) sample name i:th sample is one of the names listed by user (=if there are any TRUEs)
        if (any(seurat.objects.list[[i]]@meta.data$stim[1] == ref.samples.name.list) ){
        # if TRUE, save the i
        ref.sample.numbers <- append(ref.sample.numbers, i)
        }
    }
}else{
    ref.sample.numbers <- NULL # if no samples are listed, NULL = all pairwise anchors are found (no reference/s)
}

# 1. identify anchors using the FindIntegrationAnchors function
data.anchors <- FindIntegrationAnchors(object.list = seurat.objects.list, dims = 1:CCstocompute, anchor.features = features, reduction = anchor.identification.method, normalization.method = normalisation.method, reference = ref.sample.numbers) # dims = Which dimensions to use from the CCA to specify the neighbor search space

# 2. use these anchors to integrate the two datasets together with IntegrateData.
data.combined <- IntegrateData(anchorset = data.anchors, dims = 1:PCstocompute, normalization.method = normalisation.method) # dims = Number of PCs to use in the weighting procedure

DefaultAssay(data.combined) <- "integrated"

# Note: Skip ScaleData when using SCTransform
if (normalisation.method != "SCT"){
    data.combined <- ScaleData(data.combined, verbose = FALSE)  
}  

# Moved to next tool to clarify, there's npcs as a parameter:
# data.combined <- RunPCA(data.combined, npcs = 30, verbose = FALSE)



# Save the Robj for the next tool
save(data.combined, file="seurat_obj_combined.Robj")

## EOF



