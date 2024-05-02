# TOOL spatial-transcriptomics-seurat-spatially-variable-markvariogram-v5.R: "Seurat v5 -Identify spatially variable genes without pre-annotation" (Identify genes that have spatial patterning without taking cluster information or spatial annotation into account.)
# INPUT seurat_spatial_obj_pca.Robj: "Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL Markerplot2.pdf
# OUTPUT OPTIONAL spatially_variable_genes.tsv
# PARAMETER OPTIONAL method.to.use: "Method to use" TYPE [markvariogram: markvariogram, moransi: moransi] DEFAULT moransi (Method to use. Mark variogram takes longer to run, Morans I is faster.)
# PARAMETER OPTIONAL number.of.top.features: "Number of features to plot" TYPE INTEGER DEFAULT 6 (How many top features to plot.)
# PARAMETER OPTIONAL multi_analysis: "Multiple sample analysis" TYPE [TRUE:yes, FALSE:no] DEFAULT FALSE (If you performed multiple sample analysis i.e. if your samples have been merged/integrated into one single seurat object, select 'yes'. However, this tool will still identify spatially variable genes for each sample individually.)
# RUNTIME R-4.2.3-seurat5
# SLOTS 3
# TOOLS_BIN ""

# 2022-08-01 IH
# 2022-10-20 ML Add moransi option & output for spatially_variable_genes.tsv
# 2024-03-21 EP Update to Seurat v5

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(Biobase)

print(package.version("Seurat"))
documentVersion("Seurat", package.version("Seurat"))

# Load the R-Seurat-object (called seurat_obj)
load("seurat_spatial_obj_pca.Robj")

# Parallelise
library(future)
plan("multisession", workers = as.integer(chipster.threads.max))

if (multi_analysis) {
  # Split obejct and analyze spatially variable features separately as shown in
  # https://ucdavis-bioinformatics-training.github.io/2022-December-Spatial-Transcriptomics/data_analysis/spatial_analysis
  # Samples are analyzed seperately also because of error 'please provide the same number of observations as spatial locations'
  # in FindSpatiallyVariableFeatures() as explained by tingchiafelix in issue https://github.com/satijalab/seurat/issues/4611
  # as well as issue raised in https://github.com/satijalab/seurat/issues/8754
  experiment.slices <- SplitObject(seurat_obj, split.by = "orig.ident")
  # Remove extra images
  experiment.slices <- lapply(experiment.slices, function(sce) {
    sce@images[names(sce@images) != sce$orig.ident[1]] = NULL
    sce
  })

  # Open the pdf file for plotting
  pdf(file = "Markerplot2.pdf", width = 9, height = 12)

  for (i in 1:length(experiment.slices)) {
    seurat_obj <- FindSpatiallyVariableFeatures(experiment.slices[[i]], features = VariableFeatures(experiment.slices[[i]]), selection.method = method.to.use)
    # Visualise the identified top features
    top.features <- head(SpatiallyVariableFeatures(seurat_obj, method = method.to.use), number.of.top.features)
    print(SpatialFeaturePlot(seurat_obj, features = top.features, ncol = 3, alpha = c(0.1, 1)))
    # Print out markers into a table:
    write.table(as.matrix(top.features), file = "spatially_variable_genes.tsv", sep = "\t", row.names = T, col.names = T, quote = F)
  }

  plan("default") # return to normal

  # Close the pdf
  dev.off()


} else {
  # Find spatially variable features using markvariogram
  seurat_obj <- FindSpatiallyVariableFeatures(seurat_obj, assay = "SCT", features = VariableFeatures(seurat_obj), selection.method = method.to.use)

  plan("default") # return to normal

  # Open the pdf file for plotting
  pdf(file = "Markerplot2.pdf", width = 9, height = 12)

  # Visualise the identified top features
  top.features <- head(SpatiallyVariableFeatures(seurat_obj, method = method.to.use), number.of.top.features)

  print(SpatialFeaturePlot(seurat_obj, features = top.features, ncol = 3, alpha = c(0.1, 1)))

  # Close the pdf
  dev.off()

  # Print out markers into a table:
  write.table(as.matrix(top.features), file = "spatially_variable_genes.tsv", sep = "\t", row.names = T, col.names = T, quote = F)

}

# EOF
