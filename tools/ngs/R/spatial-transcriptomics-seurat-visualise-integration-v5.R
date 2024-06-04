# TOOL spatial-transcriptomics-seurat-visualise-integration-v5.R: "Seurat v5 -Visualize cell types after integration with scRNA-seq data" (This tool visualizes the predicted underlying composition of cell types in each spatial spot after integration with scRNA-seq reference data.)
# INPUT seurat_obj_integrated.Robj: "Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL integration_plot.pdf
# PARAMETER OPTIONAL genes: "Cell types to plot" TYPE UNCHECKED_STRING DEFAULT "L2/3 IT" (If you list multiple cell types, please use comma\(s\) \(,\) as a separator, e.g., \"L2/3 IT\,L4\".)
# PARAMETER OPTIONAL method.to.use: "Method to use" TYPE [markvariogram: markvariogram, moransi: moransi] DEFAULT moransi (Method to use. Markvariogram takes longer to run, Morans I is faster.)
# PARAMETER OPTIONAL number.of.top.features: "Number of top spatially variable cell types to plot" TYPE INTEGER DEFAULT 4 (How many top spatially variable cell types to plot.)
# PARAMETER OPTIONAL color.scale: "Determine color scale based on all cell types" TYPE [all:yes, feature:no] DEFAULT feature (Determine whether the color scale is based on all cell types or individual cell types. By default, the color scale is determined for each cell type individually and may differ between cells.)
# PARAMETER OPTIONAL multi_analysis: "Multiple sample analysis" TYPE [TRUE:yes, FALSE:no] DEFAULT FALSE (If you performed multiple sample analysis i.e. if your samples have been merged/integrated into one single seurat object, select 'yes'. This tool will visualize cell types after integration for each sample separately.)
# RUNTIME R-4.2.3-seurat5
# TOOLS_BIN ""

# 2022-08-05 IH
# 2022-10-20 ML Add moransi option

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(Biobase)

print(package.version("Seurat"))
documentVersion("Seurat", package.version("Seurat"))

# Load the R-Seurat-object (called seurat_obj)
load("seurat_obj_integrated.Robj")

DefaultAssay(seurat_obj) <- "predictions"

genes <- trimws(unlist(strsplit(genes, ",")))

if (multi_analysis) {
    # Split obejct and analyze find spatially variable features separately 
    # Samples are analyzed seperately because of error 'please provide the same number of observations as spatial locations'
    # in FindSpatiallyVariableFeatures() that is similar as explained by tingchiafelix in issue 
    # https://github.com/satijalab/seurat/issues/4611 as well as issue raised in 
    # https://github.com/satijalab/seurat/issues/8754
    experiment.slices <- SplitObject(seurat_obj, split.by = "orig.ident")

    # Remove extra images
    experiment.slices <- lapply(experiment.slices, function(sce) {
    sce@images[names(sce@images) != sce$orig.ident[1]] = NULL
    sce
    })

    # Open the pdf file for plotting
    pdf(file = "integration_plot.pdf", width = 13, height = 7)

  for (i in 1:length(experiment.slices)) {
    seurat_obj <- FindSpatiallyVariableFeatures(experiment.slices[[i]], 
        assay = "predictions", selection.method = method.to.use, 
        features = rownames(experiment.slices[[i]]), r.metric=5, slot = "data")

    top.clusters <- head(SpatiallyVariableFeatures(seurat_obj, method = method.to.use), number.of.top.features)

    # Visualise chosen features
    print(SpatialFeaturePlot(seurat_obj, features = c(genes), keep.scale = color.scale, pt.size.factor = 1.6, ncol = 2, crop = TRUE))

    # Visualise spatially variable features
    print(SpatialPlot(object = seurat_obj, features = top.clusters, keep.scale = color.scale, ncol = 2))
  }

    # Close the pdf
    dev.off()

} else {
    # Identify spatially variable features with the cell type prediction scores calculated in the integration
    seurat_obj <- FindSpatiallyVariableFeatures(seurat_obj,
        assay = "predictions", selection.method = method.to.use,
        features = rownames(seurat_obj), r.metric = 5, slot = "data", 
    )
    top.clusters <- head(SpatiallyVariableFeatures(seurat_obj, method = method.to.use), number.of.top.features)

    # Open the pdf file for plotting
    pdf(file = "integration_plot.pdf", width = 13, height = 7)

    # Visualise chosen features
    print(SpatialFeaturePlot(seurat_obj, features = c(genes), keep.scale = color.scale, pt.size.factor = 1.6, ncol = 2, crop = TRUE))
    
    # Visualise spatially variable features
    print(SpatialPlot(object = seurat_obj, features = top.clusters, keep.scale = color.scale, ncol = 2))

    # Close the pdf
    dev.off()
}

# EOF
