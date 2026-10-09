# TOOL spatial-transcriptomics-seurat-setup-v5-HD.R: "Seurat v5 -Setup and QC" (This tool sets up a Seurat object and produces quality control plots for Visium HD data. You need to give the input files as a tar package, please see the manual for details.)
# INPUT files.tar: "tar package of 10X output files" TYPE GENERIC
# OUTPUT OPTIONAL seurat_spatial_setup.Robj
# OUTPUT OPTIONAL QC_plots.pdf
# PARAMETER OPTIONAL sample_name: "Name for the sample" TYPE STRING DEFAULT "slice1" (Name for the sample. Make sure the samples are named differently if you have multiple samples.)
# PARAMETER OPTIONAL bin_sizes: "Bin sizes" TYPE STRING DEFAULT "8, 16" (List here the bin sizes, separated by comma)
# RUNTIME R-4.5.1-visium-hd
# SLOTS 4
# TOOLS_BIN ""

# 2026-02 ML 
# 2026-10-08 JV Make pdf smaller


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
library("ragg")
library("magick")

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

# There is a problem if 2um folder is included // JV
# If you comment out the --strip-components, the problem is solved // JV
system("mkdir input_folder && tar -xf files.tar -C input_folder  --strip-components=0 2> /dev/null")

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

# This pdf code makes the pdf file big > 200Mb. Underneath is png implementation that will be turned into a smaller pdf at the end.

# Open the pdf file for plotting
# pdf(file = "QC_plots.pdf", width = 13, height = 7)

# # Loop through the bins
# for (i in 1:length(assay_names)) {

#   DefaultAssay(seurat_obj) <- assay_names[i] # "Spatial.008um"
#   assay_bin <- assay_names[i]
#   nCount_bin <- paste("nCount_",assay_names[i], sep="")
#   nFeature_bin <- paste("nFeature_",assay_names[i], sep="")
#   just.bin <- sub("Spatial\\.", "", assay_names[i])

#   # vln.plot <- VlnPlot(object, features = "nCount_Spatial.008um", pt.size = 0) + theme(axis.text = element_text(size = 4)) + NoLegend()
#   vln.plot <- VlnPlot(seurat_obj, features = nCount_bin, pt.size = 0) + theme(axis.text = element_text(size = 4)) + NoLegend()
#   # count.plot <- SpatialFeaturePlot(object, features = "nCount_Spatial.008um") + theme(legend.position = "right")
#   count.plot <- SpatialFeaturePlot(seurat_obj, features = nCount_bin) + theme(legend.position = "right")

#   # note that many spots have very few counts, in-part
#   # due to low cellular density in certain tissue regions
#   print(vln.plot | count.plot)

#   seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-|^mt-|^Mt-")
#   seurat_obj[["percent.hb"]] <- PercentageFeatureSet(seurat_obj, pattern = "^Hb.*-")
#   seurat_obj[["percent.rb"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RPS|^RPL|^rps|^rpl|^Rps|^Rpl")

#   # VlnPlot(seurat_obj, features = c("nCount_Spatial", "nFeature_Spatial"), pt.size = 0.1, ncol = 2) + NoLegend()
#   # VlnPlot(seurat_obj, features = c("percent.mt", "percent.hb", "percent.rb"), pt.size = 0.1, ncol = 2) + NoLegend()
#   # SpatialFeaturePlot(seurat_obj, c("nCount_Spatial", "nFeature_Spatial", "percent.mt", "percent.rb", "percent.hb")) # + theme(legend.position = "right")
#   print(VlnPlot(seurat_obj, features = c(nCount_bin, nFeature_bin), pt.size = 0.1, ncol = 2) + NoLegend())
#   print(VlnPlot(seurat_obj, features = c("percent.mt", "percent.hb"), pt.size = 0.1, ncol = 2) + NoLegend()) # , "percent.rb"
#   print(SpatialFeaturePlot(seurat_obj, c(nCount_bin, nFeature_bin, "percent.mt", "percent.hb"))) # + theme(legend.position = "right") , "percent.rb"


#   # Top expressing genes
#   # Code from: https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/seurat/seurat_07_spatial.html#Top_expressed_genes)
#   # C <- LayerData(seurat_obj, assay = "Spatial.008um", layer = "counts")
#   C <- LayerData(seurat_obj, assay = assay_bin, layer = "counts")
#   C@x <- C@x / rep.int(colSums(C), diff(C@p))
#   most_expressed <- order(Matrix::rowSums(C), decreasing = T)[20:1]
#   print(boxplot(as.matrix(t(C[most_expressed, ])),
#     cex = 0.1, las = 1, xlab = "% total count per spot",
#     col = (scales::hue_pal())(20)[20:1], horizontal = TRUE
#   ))
# } # end looping for bins = assays

# # Close the pdf
# dev.off()


# Plot png and turn them into a pdf at the end.


output_dir <- paste0(getwd(), "/QC_plots")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)


for (i in seq_along(assay_names)) {

  assay_bin <- assay_names[i]
  DefaultAssay(seurat_obj) <- assay_bin
  
  nCount_bin <- paste0("nCount_", assay_bin)
  nFeature_bin <- paste0("nFeature_", assay_bin)
  just.bin <- sub("Spatial\\.", "", assay_bin)

  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-|^mt-|^Mt-", assay = assay_bin)
  seurat_obj[["percent.hb"]] <- PercentageFeatureSet(seurat_obj, pattern = "^HB[^(P)]|^Hb[^(p)]", assay = assay_bin)
  seurat_obj[["percent.rb"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RPS|^RPL|^rps|^rpl|^Rps|^Rpl", assay = assay_bin)


  vln.plot <- VlnPlot(seurat_obj, features = nCount_bin, pt.size = 0) + theme(axis.text = element_text(size = 4)) + NoLegend()
  count.plot <- SpatialFeaturePlot(seurat_obj, features = nCount_bin) + theme(legend.position = "right")

  p0 <- print(vln.plot | count.plot)

  # Create the four plots
  p1 <- (VlnPlot(seurat_obj, features = c(nCount_bin, nFeature_bin), pt.size = 0.1, ncol = 2) + NoLegend())

  
  p2 <- (VlnPlot(seurat_obj, features = c("percent.mt", "percent.hb"), pt.size = 0.1, ncol = 2) + NoLegend()) # , "percent.rb"

  
  p3 <- (SpatialFeaturePlot(seurat_obj, c(nCount_bin, nFeature_bin, "percent.mt", "percent.hb"))) # + theme(legend.position = "right") , "percent.rb"

  C <- LayerData(seurat_obj, assay = assay_bin, layer = "counts")
  C@x <- C@x / rep.int(colSums(C), diff(C@p))
  most_expressed <- order(Matrix::rowSums(C), decreasing = T)[20:1]
  top_genes <- rownames(C)[most_expressed]


#  p4 <- (boxplot(as.matrix(t(C[most_expressed, ])),
#      cex = 0.1, las = 1, xlab = "% total count per spot",
#      col = (scales::hue_pal())(20)[20:1], horizontal = TRUE
#    ))
  # ggplot version of p4, this needs to be rechecked, so far looks good.

  top_df <- data.frame(
    gene = factor(rep(top_genes, each = ncol(C)), levels = top_genes),
    pct  = as.vector(t(as.matrix(C[most_expressed, ])))
  )
  rm(C)
  
  p4 <- ggplot(top_df, aes(x = pct, y = gene, fill = gene)) +
  geom_boxplot(outlier.size = 0.1) +
  scale_fill_manual(values = rev(scales::hue_pal()(20))) +
  labs(x = "% total count per bin", y = NULL,
    title = paste("Top 20 expressed genes,", just.bin)) +
  theme_classic() +
  NoLegend()
  
  # Put plots into a list
  plots <- list(p0, p1, p2, p3, p4)
  
  # Loop over the plots and save each as a separate PNG
  for (j in seq_along(plots)) {
    
    filename <- file.path(
      output_dir,
      paste0(just.bin, "_plot_", j, ".png")
    )
    
    agg_png(
      filename = filename,
      width = 10,
      height = 10,
      res = 300,
      bitsize = 16,
      units = "in"
    )
    
    print(plots[[j]])
    
    dev.off()
  }
}

png_files <- list.files(
  output_dir,
  pattern = "\\.png$",
  full.names = TRUE
)

# Sort them so the PDF follows assay -> plot order
png_files <- sort(png_files)

# Read all PNGs
images <- image_read(png_files)

image_write(images, format="pdf", path = "QC_plots.pdf")

# Save the Robj for the next tool
save(seurat_obj, file = "seurat_spatial_setup.Robj")


# EOF