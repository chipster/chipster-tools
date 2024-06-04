# TOOL spatial-transcriptomics-seurat-spatially-variable-markvariogram-v5.R: "Seurat v5 -Identify spatially variable genes without pre-annotation" (This tool identifies spatially variable genes without cluster annotation information and visualizes these genes on top of the tissue image.)
# INPUT seurat_spatial_obj_pca.Robj: "Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL Markerplot2.pdf
# OUTPUT OPTIONAL spatially_variable_genes.tsv
# PARAMETER OPTIONAL method.to.use: "Method to use" TYPE [markvariogram: markvariogram, moransi: moransi] DEFAULT moransi (Method to use. Mark variogram takes longer to run, Moran's I is faster.)
# PARAMETER OPTIONAL number_of_var_feats: "Number of variable genes to use" TYPE INTEGER DEFAULT 1000 (Number of variable genes to use for identifying highest spatially variable genes. You can speed up the computation by choosing a smaller number of variable genes. This number should be less than or equal to the number of variable genes in the Seurat object.)
# PARAMETER OPTIONAL number.of.top.features: "Number of genes to plot" TYPE INTEGER DEFAULT 6 (How many top features to plot.)
# PARAMETER OPTIONAL color.scale: "Determine color scale based on all genes" TYPE [all:yes, feature:no] DEFAULT feature (Determine whether the color scale is based on all genes or individual genes. By default, the color scale is determined for each gene individually and may differ between genes.)
# PARAMETER OPTIONAL multi_analysis: "Multiple sample analysis" TYPE [TRUE:yes, FALSE:no] DEFAULT FALSE (If you performed multiple sample analysis i.e. if your samples have been merged/integrated into one single seurat object, select 'yes'. This tool will identify spatially variable genes for each sample individually.)
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

  # Print warnings after SplitObject() but suppress them from Chipster note
  experiment.slices <- tryCatch({
    SplitObject(seurat_obj, split.by = "orig.ident")
  }, warning=function(w) {
    print(w) 
    return(suppressWarnings(SplitObject(seurat_obj, split.by = "orig.ident")))
  })

  # Remove extra images
  experiment.slices <- lapply(experiment.slices, function(sce) {
    sce@images[names(sce@images) != sce$orig.ident[1]] = NULL
    sce
  })

  # Open the pdf file for plotting
  pdf(file = "Markerplot2.pdf", width = 9, height = 12)

  for (i in 1:length(experiment.slices)) {
    # Number of variable features in Seurat object
    variable_feats <- VariableFeatures(experiment.slices[[i]])

    # Use number of selected variable features 
    if (number_of_var_feats > length(variable_feats)) {
      stop(paste("CHIPSTER-NOTE: ", "You have selected", number_of_var_feats, "to use as the number of variable features, but there are only", length(variable_feats), "in the object. Please select less than or equal to", length(variable_feats),"variable genes."))
    }
    variable_feats <- variable_feats[1:number_of_var_feats]
    
    seurat_obj <- FindSpatiallyVariableFeatures(experiment.slices[[i]], assay = "SCT", features = variable_feats, selection.method = method.to.use)

    # Visualise the identified top features
    top.features <- head(SpatiallyVariableFeatures(seurat_obj, method = method.to.use), number.of.top.features)
    print(SpatialFeaturePlot(seurat_obj, features = top.features, keep.scale = color.scale, ncol = 3, alpha = c(0.1, 1)))
    # Print out markers into a table:
    write.table(as.matrix(top.features), file = "spatially_variable_genes.tsv", sep = "\t", row.names = T, col.names = T, quote = F)
  }

  # Return to normal
  plan("default") 

  # Close the pdf
  dev.off()

} else {
  # Number of variable features in Seurat object
  variable_feats <- VariableFeatures(seurat_obj)

  # Use number of selected variable features 
  if (number_of_var_feats > length(variable_feats)) {
    stop(paste("CHIPSTER-NOTE: ", "You have selected", number_of_var_feats, "to use as the number of variable features, but there are only", length(variable_feats), "in the object. Please select less than or equal to", length(variable_feats),"variable genes."))
  }
  variable_feats <- variable_feats[1:number_of_var_feats]

  seurat_obj <- FindSpatiallyVariableFeatures(seurat_obj, assay = "SCT", features = variable_feats, selection.method = method.to.use)

  plan("default") # return to normal

  # Open the pdf file for plotting
  pdf(file = "Markerplot2.pdf", width = 9, height = 12)

  # Visualise the identified top features
  top.features <- head(SpatiallyVariableFeatures(seurat_obj, method = method.to.use), number.of.top.features)

  print(SpatialFeaturePlot(seurat_obj, features = top.features, keep.scale = color.scale, ncol = 3, alpha = c(0.1, 1)))

  # Close the pdf
  dev.off()

  # Print out markers into a table:
  write.table(as.matrix(top.features), file = "spatially_variable_genes.tsv", sep = "\t", row.names = T, col.names = T, quote = F)

}

# EOF
