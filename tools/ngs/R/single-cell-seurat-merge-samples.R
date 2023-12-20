# TOOL single-cell-seurat-merge-samples.R: "Seurat v5 -Merge samples" (Merge multiple samples for joined analysis.)
# INPUT samples{...}.Robj: "Samples to combine" TYPE GENERIC
# OUTPUT seurat_obj_combined.Robj
# OUTPUT OPTIONAL Dispersion_plot.pdf
# OUTPUT OPTIONAL PCAplots.pdf
# OUTPUT OPTIONAL PCAloadings.txt
# PARAMETER OPTIONAL normalisation.method: "Normalization method to perform" TYPE [LogNormalize:"Global scaling normalization", SCT:"SCTransform"] DEFAULT LogNormalize (Normalize data with global scaling normalization or SCTransform.)
# PARAMETER OPTIONAL totalexpr: "Scaling factor in the normalization" TYPE INTEGER DEFAULT 10000 (Scale each cell to this total number of transcripts.)
# PARAMETER OPTIONAL num.features: "Number of variable genes to return" TYPE INTEGER DEFAULT 2000 (Number of features to select as top variable features, i.e. how many features returned.)
# PARAMETER OPTIONAL num.of.pcas: "Number of PCs to compute" TYPE INTEGER DEFAULT 30 (How many principal components to compute and store. If you get an error message, try lowering the number. This might happen especially if you have low number of cells in your data.)
# PARAMETER OPTIONAL num.of.heatmaps: "Number of principal components to plot as heatmaps" TYPE INTEGER DEFAULT 12 (How many principal components to plot as heatmaps.)
# PARAMETER OPTIONAL loadings: "Print loadings in a file" TYPE [TRUE: yes, FALSE: no] DEFAULT FALSE (Print the PC loadings to a txt file.)
# PARAMETER OPTIONAL num.of.genes.loadings: "Number of genes to list in the loadings file" TYPE INTEGER DEFAULT 5 (How many genes to list in the loadings txt file.)
# RUNTIME R-4.3.2-single-cell
# SLOTS 3
# TOOLS_BIN ""

# 2023-11-29 IH

library(Seurat)
library(dplyr)
library(Matrix)
library(gplots)
library(ggplot2)
options(Seurat.object.assay.version = "v5")

# Load the R-Seurat-objects (called seurat_obj -that's why we need to rename them here)
input.names <- read.table("chipster-inputs.tsv", header = FALSE, sep = "\t")
for (i in 1:nrow(input.names)) {
    load(input.names[i, 1])
    name.of.obj <- paste("seurat_obj_", i, sep = "")
    assign(name.of.obj, seurat_obj)
}

seurat.objects.list <- as.list(mget(objects(pattern = "seurat_obj_")))

#if (normalisation.method == "SCT") {
    #features <- SelectIntegrationFeatures(seurat.objects.list)
    #seurat.objects.list <- lapply(X = seurat.objects.list, FUN = SCTransform)
#}
# merge first seurat object in the list with the other seurat object(s)
seurat_obj <- merge(x = seurat.objects.list[[1]], y = seurat.objects.list[-1])
seurat_obj

# Normalisation
if (normalisation.method == "SCT") {
    seurat_obj <- SCTransform(seurat_obj, assay = "RNA", new.assay.name = "SCT", vars.to.regress = c("percent.mt"), variable.features.n = num.features, verbose = FALSE)
    #VariableFeatures(seurat_obj) <- features
    #seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = num.features)
    #VariableFeatures(seurat_obj) <- features 
} else if (normalisation.method == "LogNormalize") {
    seurat_obj <- NormalizeData(object = seurat_obj, normalization.method = "LogNormalize", scale.factor = totalexpr)
    # Detection of variable genes across the single cells
    seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = num.features)

    ## Scaling:
    seurat_obj <- ScaleData(object = seurat_obj, vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = FALSE)
    print(seurat_obj)
}

# Open the pdf file for plotting
pdf(file = "Dispersion_plot.pdf", width = 13, height = 7)

# Dispersion plot:
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seurat_obj), 10)
# Plot variable features with and without labels
dim(VariableFeatures(seurat_obj))
plot1 <- VariableFeaturePlot(seurat_obj)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))

textplot(paste("\v \v Number of \n \v \v variable \n \v \v genes: \n \v \v", length(VariableFeatures(object = seurat_obj)), " \n  \n \v \v Number of \n \v \v cells: \n \v \v", length(colnames(x = seurat_obj))), halign = "center", valign = "center", cex = 2)
dev.off() # close the pdf


# Save the Robj for the next tool
save(seurat_obj, file = "seurat_obj_combined.Robj")

## EOF
