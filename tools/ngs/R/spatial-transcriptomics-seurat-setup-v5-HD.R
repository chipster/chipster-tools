# TOOL spatial-transcriptomics-seurat-setup-v5-HD.R: "Seurat v5 HD -Setup and QC" (This tool sets up a Seurat object for spatial transcriptomics data.)
# INPUT files.tar: "tar package of 10X output files" TYPE GENERIC
# OUTPUT OPTIONAL seurat_spatial_setup.Robj
# OUTPUT OPTIONAL QC_plots.pdf
# PARAMETER OPTIONAL sample_name: "Name for the sample" TYPE STRING DEFAULT "slice1" (Name for the sample. Make sure the samples are named differently if you have multiple samples.)
# PARAMETER OPTIONAL bin_sizes: "Bin sizes" TYPE STRING DEFAULT "8, 16" (List here the bin sizes, separated by comma)
# RUNTIME R-4.5.1-seurat5
# SLOTS 3
# TOOLS_BIN ""

# 2026-02 ML 

# RUNTIME R-4.2.3-seurat5 -> R-4-5-1-seurat5 

# install.packages('ggplot2', repos='http://cran.us.r-project.org') # 3.5.2
# install.packages('patchwork', repos='http://cran.us.r-project.org') # 1.3.1 
# install.packages('hdf5r', repos='http://cran.us.r-project.org') # 1.3.10 / 1.3.12
# install.packages('arrow', repos='http://cran.us.r-project.org') # 20.0.0.2
# install.packages('Seurat', repos='http://cran.us.r-project.org') # 5.3.0

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(Matrix)
library(Biobase)

# Get the Seurat version info for displaying:
source(file.path(chipster.common.lib.path, "tool-utils.R"))
print(package.version("Seurat"))
documentVersion("Seurat", package.version("Seurat"))


# Note for developers:
# More elegant file handling in non-HD version. Might not work here, as there are identically named files and subfolders 
# within the different "bin" folders. 


# Open tar package into a new folder called input_folder. 
# Keep the directory structure!
# Note: now we are assuming the binned_outputs folder is within a user defined folder.
# -C changes to the specified directory before unpacking. 
# --strip-components 1 removes 1 directory from the filenames stored in the archive. 
# "2> /dev/null" can be used to redirect the errors to /dev/null (Mac OS X uses BSD tar and creates some extra info that is not recognized by GNU tar which causes messages: tar: Ignoring unknown extended header keyword 'SCHILY.fflags')
system("mkdir input_folder && tar -xf files.tar -C input_folder --strip-components=2 2> /dev/null")

# For testing:
# die here:
# charToRaw(testi2)

# Replace empty spaces in sample name with underscore _ :
sample_name <- gsub(" ", "_", sample_name)

# bin sizes:
# bin_sizes <-  "8, 16"
bin_sizes <- as.numeric(trimws(unlist(strsplit(bin_sizes, ","))))
print(bin_sizes)

# Load spatial data, returns a Seurat object
# slice = name for the stored image of the tissue slice later used in the analysis
seurat_obj <- Load10X_Spatial(data.dir = "input_folder/", assay = "Spatial", slice = sample_name, bin.size = bin_sizes) # bin.size = c(8, 16)

# # Sometimes the coordinates of the spatial image are characters instead of integers.
# # Check this and switch them to intergers if need be:
# if (class(seurat_obj@images[[sample_name]]@coordinates[["tissue"]]) == "character") {
#   seurat_obj@images[[sample_name]]@coordinates[["tissue"]] <- as.integer(seurat_obj@images[[sample_name]]@coordinates[["tissue"]])
#   seurat_obj@images[[sample_name]]@coordinates[["row"]] <- as.integer(seurat_obj@images[[sample_name]]@coordinates[["row"]])
#   seurat_obj@images[[sample_name]]@coordinates[["col"]] <- as.integer(seurat_obj@images[[sample_name]]@coordinates[["col"]])
#   seurat_obj@images[[sample_name]]@coordinates[["imagerow"]] <- as.integer(seurat_obj@images[[sample_name]]@coordinates[["imagerow"]])
#   seurat_obj@images[[sample_name]]@coordinates[["imagecol"]] <- as.integer(seurat_obj@images[[sample_name]]@coordinates[["imagecol"]])
# }

# Add samplename to metadata field (for plotting)
# seurat_obj <- RenameIdents(seurat_obj, "SeuratProject" = sample_name)
seurat_obj@meta.data$orig.ident <- sample_name


assay_names <- Assays(seurat_obj) # [1] "Spatial.008um" "Spatial.016um"

# Open the pdf file for plotting
pdf(file = "QC_plots.pdf", width = 13, height = 7)

# Loop through the bins
for (i in 1:length(assay_names)) {

  DefaultAssay(seurat_obj) <- assay_names[i] # "Spatial.008um"
  assay_bin <- assay_names[i]
  nCount_bin <- paste("nCount_",assay_names[i], sep="")
  nFeature_bin <- paste("nFeature_",assay_names[i], sep="")
  just.bin <- sub("Spatial\\.", "", assay_names[i])

  # vln.plot <- VlnPlot(object, features = "nCount_Spatial.008um", pt.size = 0) + theme(axis.text = element_text(size = 4)) + NoLegend()
  vln.plot <- VlnPlot(seurat_obj, features = nCount_bin, pt.size = 0) + theme(axis.text = element_text(size = 4)) + NoLegend()
  # count.plot <- SpatialFeaturePlot(object, features = "nCount_Spatial.008um") + theme(legend.position = "right")
  count.plot <- SpatialFeaturePlot(seurat_obj, features = nCount_bin) + theme(legend.position = "right")

  # note that many spots have very few counts, in-part
  # due to low cellular density in certain tissue regions
  print(vln.plot | count.plot)

  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-|^mt-|^Mt-")
  seurat_obj[["percent.hb"]] <- PercentageFeatureSet(seurat_obj, pattern = "^Hb.*-")
  seurat_obj[["percent.rb"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RPS|^RPL|^rps|^rpl|^Rps|^Rpl")

  # VlnPlot(seurat_obj, features = c("nCount_Spatial", "nFeature_Spatial"), pt.size = 0.1, ncol = 2) + NoLegend()
  # VlnPlot(seurat_obj, features = c("percent.mt", "percent.hb", "percent.rb"), pt.size = 0.1, ncol = 2) + NoLegend()
  # SpatialFeaturePlot(seurat_obj, c("nCount_Spatial", "nFeature_Spatial", "percent.mt", "percent.rb", "percent.hb")) # + theme(legend.position = "right")
  print(VlnPlot(seurat_obj, features = c(nCount_bin, nFeature_bin), pt.size = 0.1, ncol = 2) + NoLegend())
  print(VlnPlot(seurat_obj, features = c("percent.mt", "percent.hb"), pt.size = 0.1, ncol = 2) + NoLegend()) # , "percent.rb"
  print(SpatialFeaturePlot(seurat_obj, c(nCount_bin, nFeature_bin, "percent.mt", "percent.hb"))) # + theme(legend.position = "right") , "percent.rb"


  # Top expressing genes
  # Code from: https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/seurat/seurat_07_spatial.html#Top_expressed_genes)
  # C <- LayerData(seurat_obj, assay = "Spatial.008um", layer = "counts")
  C <- LayerData(seurat_obj, assay = assay_bin, layer = "counts")
  C@x <- C@x / rep.int(colSums(C), diff(C@p))
  most_expressed <- order(Matrix::rowSums(C), decreasing = T)[20:1]
  print(boxplot(as.matrix(t(C[most_expressed, ])),
    cex = 0.1, las = 1, xlab = "% total count per spot",
    col = (scales::hue_pal())(20)[20:1], horizontal = TRUE
  ))
} # end looping for bins = assays

# Close the pdf
dev.off()

# Save the Robj for the next tool
save(seurat_obj, file = "seurat_spatial_setup.Robj")

# EOF
