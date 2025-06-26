# TOOL spatial-transcriptomics-seurat-setup-v5.R: "Seurat v5 -Setup and QC" (This tool sets up a Seurat object for spatial transcriptomics data.)
# INPUT OPTIONAL files.tar: "tar package of 10X output files" TYPE GENERIC
# OUTPUT OPTIONAL seurat_spatial_setup.Robj
# OUTPUT OPTIONAL QC_plots.pdf
# PARAMETER OPTIONAL sample_name: "Name for the sample" TYPE STRING DEFAULT "slice1" (Name for the sample. Make sure the samples are named differently if you have multiple samples.)
# RUNTIME R-4.2.3-seurat5
# SLOTS 2
# TOOLS_BIN ""

# 2022-07-15 IH
# 2022-10-13 ML Coordinates to integers -check and input folder handling
# 2022-10-18 ML Add samplename to a metadata field
# 2022-11-01 ML Add top expressed genes boxplot
# 2023-09-08 IH add percent.rb to QC
# 2024-03-21 EP Update to Seurat v5
# 2025-06-26 ML Add HB plotting on top of tissue image

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(Matrix)
library(Biobase)

source(file.path(chipster.common.lib.path, "tool-utils.R"))
print(package.version("Seurat"))
documentVersion("Seurat", package.version("Seurat"))

# Replace empty spaces in sample name with underscore _ :
sample_name <- gsub(" ", "_", sample_name)

# Read the contents of the tar file into a list
system("tar tf files.tar > tar.contents 2>> log.txt")
file.list <- scan("tar.contents", what = "", sep = "\n")

# Check that the input is a valid tar file
if (length(file.list) == 0) {
  stop(paste("CHIPSTER-NOTE: ", "It seems your input file is not a valid Tar package. Please check your input file."))
}

# Make an output folder
system("mkdir output_folder")
system("mkdir output_folder/spatial")

# Open tar package. Make a folder called input_folder, open the tar there so that each file
# will be on the root level (remove everything from the name until the last "/" with the --xform option)

# system("mkdir input_folder; cd input_folder; tar xf ../files.tar")
system("mkdir input_folder; cd input_folder; tar xf ../files.tar --xform='s#^.+/##x' 2>> log.txt")

# Code for manipulating aggregation files (not needed atm)
# if (file.exists("aggregation.csv") & file.exists("aggr_tissue_positions_list.csv")) {
# aggregation splitting step based on the aggregation.csv
# aggregation.csv, aggr_tissue_position_list.csv -> tissue_positions_list.csv }

# Rename and move filtered_feature_bc_matrix.h5
file1 <- list.files("input_folder", pattern = "feature_bc_matrix.h5", full.names = TRUE)
file1 <- paste("mv", file1, "output_folder/filtered_feature_bc_matrix.h5")
system(file1)

# Rename and move tissue_positions_list.csv
tissue_positions <- list.files("input_folder", pattern = "tissue_positions", full.names = TRUE)
if (file.exists(tissue_positions)) {
  tissue_positions <- paste("mv", tissue_positions, "output_folder/spatial/tissue_positions_list.csv")
  system(tissue_positions)
}

# Copy image files to spatial subdir
file.copy("input_folder/scalefactors_json.json", "output_folder/spatial")
file.copy("input_folder/tissue_hires_image.png", "output_folder/spatial")
file.copy("input_folder/tissue_lowres_image.png", "output_folder/spatial")
file.copy("input_folder/tissue_positions.csv", "output_folder/spatial/tissue_positions_list.csv")

# Load spatial data, returns a Seurat object
# slice = name for the stored image of the tissue slice later used in the analysis
seurat_obj <- Load10X_Spatial(data.dir = "output_folder", assay = "Spatial", slice = sample_name)

# Sometimes the coordinates of the spatial image are characters instead of integers.
# Check this and switch them to intergers if need be:
if (class(seurat_obj@images[[sample_name]]@coordinates[["tissue"]]) == "character") {
  seurat_obj@images[[sample_name]]@coordinates[["tissue"]] <- as.integer(seurat_obj@images[[sample_name]]@coordinates[["tissue"]])
  seurat_obj@images[[sample_name]]@coordinates[["row"]] <- as.integer(seurat_obj@images[[sample_name]]@coordinates[["row"]])
  seurat_obj@images[[sample_name]]@coordinates[["col"]] <- as.integer(seurat_obj@images[[sample_name]]@coordinates[["col"]])
  seurat_obj@images[[sample_name]]@coordinates[["imagerow"]] <- as.integer(seurat_obj@images[[sample_name]]@coordinates[["imagerow"]])
  seurat_obj@images[[sample_name]]@coordinates[["imagecol"]] <- as.integer(seurat_obj@images[[sample_name]]@coordinates[["imagecol"]])
}

# Add samplename to metadata field (for plotting)
seurat_obj <- RenameIdents(seurat_obj, "SeuratProject" = sample_name)
seurat_obj@meta.data$orig.ident <- sample_name

# Open the pdf file for plotting
pdf(file = "QC_plots.pdf", width = 13, height = 7)

seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-|^mt-|^Mt-")
seurat_obj[["percent.hb"]] <- PercentageFeatureSet(seurat_obj, pattern = "^Hb.*-")
seurat_obj[["percent.rb"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RPS|^RPL|^rps|^rpl|^Rps|^Rpl")

VlnPlot(seurat_obj, features = c("nCount_Spatial", "nFeature_Spatial"), pt.size = 0.1, ncol = 2) + NoLegend()
VlnPlot(seurat_obj, features = c("percent.mt", "percent.hb", "percent.rb"), pt.size = 0.1, ncol = 2) + NoLegend()
SpatialFeaturePlot(seurat_obj, c("nCount_Spatial", "nFeature_Spatial", "percent.mt", "percent.rb", "percent.hb")) # + theme(legend.position = "right")

# Top expressing genes
# Code from: https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/seurat/seurat_07_spatial.html#Top_expressed_genes)
C <- LayerData(seurat_obj, assay = "Spatial", layer = "counts")
C@x <- C@x / rep.int(colSums(C), diff(C@p))
most_expressed <- order(Matrix::rowSums(C), decreasing = T)[20:1]
boxplot(as.matrix(t(C[most_expressed, ])),
  cex = 0.1, las = 1, xlab = "% total count per spot",
  col = (scales::hue_pal())(20)[20:1], horizontal = TRUE
)

# Close the pdf
dev.off()

# Save the Robj for the next tool
save(seurat_obj, file = "seurat_spatial_setup.Robj")

# EOF
