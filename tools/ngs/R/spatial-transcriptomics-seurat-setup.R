# TOOL spatial-transcriptomics-seurat-setup.R: "Seurat v4 -Setup and QC" (Setup the Seurat object from spatial data.)
# INPUT OPTIONAL files.tar: "tar package of 10X output files" TYPE GENERIC
# OUTPUT OPTIONAL seurat_spatial_setup.Robj
# OUTPUT OPTIONAL QC_plots.pdf 
# OUTPUT OPTIONAL log.txt
# PARAMETER OPTIONAL sample_name: "Sample name" TYPE STRING DEFAULT Sample_A (Name the sample.)
# RUNTIME R-4.1.0-single-cell

library(Seurat)
library(SeuratData)
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
system("mkdir input_folder; cd input_folder; tar xf ../files.tar --xform='s#^.+/##x' 2>> log.txt")
list.files("input_folder")

if (file.exists("aggregation.csv") & file.exists("aggr_tissue_positions_list.csv")) {
    #aggregation splitting step based on the aggregation.csv with bishwa's code
    #aggregation.csv, aggr_tissue_position_list.csv -> tissue_positions_list.csv
    #file.copy("tissue_positions_list.csv", "output_folder/spatial")
    #counter=1
    #for i in {Sample_A,Sample_B}
    #do
	  #echo -e $i"\t"$counter
	  #grep "\-"$counter aggr_tissue_positions_list.csv > spatial/$i/spatial/tissue_positions_list.csv
	  #counter=$((counter+1))
    #done
} 

system("mv input_folder/Visium_Adult_Mouse_Brain_filtered_feature_bc_matrix.h5 input_folder/filtered_feature_bc_matrix.h5") #for now

# copy image files to spatial subdir
file.copy("input_folder/scalefactors_json.json", "output_folder/spatial")
file.copy("input_folder/tissue_hires_image.png", "output_folder/spatial")
file.copy("input_folder/tissue_lowres_image.png", "output_folder/spatial")
file.copy("input_folder/tissue_positions_list.csv", "output_folder/spatial") #before aggregation works
file.copy("input_folder/filtered_feature_bc_matrix.h5", "output_folder")

#Load spatial data, returns a Seurat object
#slice = name for the stored image of the tissue slice later used in the analysis
#if multiple samples needs to be named uniquely
seurat_obj <- Load10X_Spatial(data.dir = "output_folder", assay = "Spatial", slice = "slice1")

# Open the pdf file for plotting
pdf(file="QC_plots.pdf", width=13, height=7) 

 plot1 <- VlnPlot(seurat_obj, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
 plot2 <- SpatialFeaturePlot(seurat_obj, features = "nCount_Spatial") + theme(legend.position = "right")
 CombinePlots(plots = list(plot1, plot2))

# close the pdf
 dev.off() 


# Save the Robj for the next tool
save(seurat_obj, file = "seurat_spatial_setup.Robj")

#EOF
