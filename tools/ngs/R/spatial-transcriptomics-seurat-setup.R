# TOOL spatial-transcriptomics-seurat-setup.R: "Seurat v4 -Setup and QC" (Setup the Seurat object from spatial data.)
# INPUT OPTIONAL files.tar: "tar package of 10X output files" TYPE GENERIC
# OUTPUT OPTIONAL seurat_spatial_setup.Robj
# OUTPUT OPTIONAL QC_plots.pdf 
# PARAMETER OPTIONAL sample_name: "Name for the sample" TYPE STRING DEFAULT "slice1" (Name for the sample in multisample analysis.)
# RUNTIME R-4.1.0-single-cell
# SLOTS 2

# 2022-07-15 IH

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(Matrix)

# Read the contents of the tar file into a list
system("tar tf files.tar > tar.contents 2>> log.txt")
file.list <- scan("tar.contents", what = "", sep = "\n")

# Check that the input is a valid tar file
  if (length(file.list) == 0) {
    stop(paste("CHIPSTER-NOTE: ", "It seems your input file is not a valid Tar package. Please check your input file."))
  }

# make an output folder
system("mkdir output_folder")
system("mkdir output_folder/spatial")

# Open tar package. Make a folder called input_folder, open the tar there so that each file
# will be on the root level (remove everything from the name until the last "/" with the --xform option)
system("mkdir input_folder; cd input_folder; tar xf ../files.tar")

# code for manipulating aggregation files (not needer atm)
#if (file.exists("aggregation.csv") & file.exists("aggr_tissue_positions_list.csv")) {
    #aggregation splitting step based on the aggregation.csv
    #aggregation.csv, aggr_tissue_position_list.csv -> tissue_positions_list.csv }

#rename and move filtered_feature_bc_matrix.h5
file1 <- list.files("input_folder", pattern = "feature_bc_matrix.h5", full.names = TRUE)
file1 <- paste("mv",file1,"output_folder/filtered_feature_bc_matrix.h5") 
system(file1)

#rename and move tissue_positions_list.csv 
tissue_positions <- list.files("input_folder", pattern = "tissue_positions", full.names = TRUE)
if (file.exists(tissue_positions)) {
        tissue_positions <- paste("mv",tissue_positions,"output_folder/spatial/tissue_positions_list.csv") 
        system(tissue_positions)
}

# copy image files to spatial subdir
file.copy("input_folder/scalefactors_json.json", "output_folder/spatial")
file.copy("input_folder/tissue_hires_image.png", "output_folder/spatial")
file.copy("input_folder/tissue_lowres_image.png", "output_folder/spatial")
file.copy("input_folder/tissue_positions.csv", "output_folder/spatial/tissue_positions_list.csv")

#Load spatial data, returns a Seurat object
#slice = name for the stored image of the tissue slice later used in the analysis
seurat_obj <- Load10X_Spatial(data.dir = "output_folder", assay = "Spatial", slice = sample_name)


# Open the pdf file for plotting
pdf(file="QC_plots.pdf", width=13, height=7)

seurat_obj <- PercentageFeatureSet(seurat_obj, "^mt-", col.name = "percent_mito")
seurat_obj <- PercentageFeatureSet(seurat_obj, "^Hb.*-", col.name = "percent_hb")

VlnPlot(seurat_obj, features = c("nCount_Spatial",  "nFeature_Spatial"), pt.size = 0.1, ncol = 2) + NoLegend()
VlnPlot(seurat_obj, features = c("percent_mito","percent_hb"), pt.size = 0.1, ncol = 2) + NoLegend()
SpatialFeaturePlot(seurat_obj, c("nCount_Spatial", "nFeature_Spatial", "percent_mito")) + theme(legend.position = "right")

# close the pdf
 dev.off() 

# Save the Robj for the next tool
save(seurat_obj, file = "seurat_spatial_setup.Robj")

#EOF