# TOOL spatial-transcriptomics-seurat-QC.R: "Seurat v4 -Setup and QC" (Setup the Seurat object, make quality control plots and filter out genes.)
# INPUT OPTIONAL files.tar: "tar package of 10X output files" TYPE GENERIC
# OUTPUT OPTIONAL QC_plot.pdf 
# OUTPUT OPTIONAL seurat_spatial_setup.Robj
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

#open tar into datadir
system("mkdir datadir; cd datadir; tar xf ../files.tar --xform='s#^.+/##x' 2>> log.txt")

#Load spatial data, returns a Seurat object
image.file <- Read10X_Image("datadir", image.name = "tissue_lowres_image.png")
seurat_obj <- Load10X_Spatial(data.dir = "datadir", filename = "Visium_Adult_Mouse_Brain_filtered_feature_bc_matrix.h5", assay = "Spatial", image = image.file)

# Open the pdf file for plotting
pdf(file="QC_plot.pdf", width=13, height=7) 

 plot1 <- VlnPlot(seurat_obj, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
 plot2 <- SpatialFeaturePlot(seurat_obj, features = "nCount_Spatial") + theme(legend.position = "right")
 CombinePlots(plots = list(plot1, plot2))

# # close the pdf
 dev.off() 

# Save the Robj for the next tool
save(seurat_obj, file="seurat_spatial_setup.Robj")


## EOF