# TOOL spatial-transcriptomics-seurat-subset-cells-v5.R: "Seurat v5 -Subset out cells" (Subset out anatomical regions based on exact positions.)
# INPUT seurat_obj_subset1.Robj: "Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL seurat_obj_subset2.Robj
# OUTPUT OPTIONAL subset2.pdf
# PARAMETER OPTIONAL rows: "Choose cell rows" TYPE INTEGER DEFAULT 400 (Check with spatialdimplot which cell rows to remove.)
# PARAMETER OPTIONAL columns: "Choose cell columns" TYPE INTEGER DEFAULT 150 (Check with spatialdimplot which cell columns to remove.)
# RUNTIME R-4.2.3-single-cell
# TOOLS_BIN ""

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

# Load the R-Seurat-object (called seurat_obj)
load("seurat_obj_subset1.Robj")

sample_name <- Images(seurat_obj)

seurat_obj@images[[sample_name]]@coordinates <- subset(seurat_obj@images[[sample_name]]@coordinates, imagerow > rows | imagecol < columns, invert = TRUE)
# cortex <- subset(cortex, anterior1_imagerow > 275 & anterior1_imagecol > 370, invert = TRUE)
# cortex <- subset(cortex, anterior1_imagerow > 250 & anterior1_imagecol > 440, invert = TRUE)

# Open the pdf file for plotting
pdf(file = "subset2.pdf", width = 13, height = 7)

p1 <- SpatialDimPlot(seurat_obj, crop = TRUE, label = TRUE)
p2 <- SpatialDimPlot(seurat_obj, crop = FALSE, label = TRUE, pt.size.factor = 1, label.size = 3)
p1 + p2

# close the pdf
dev.off()

# Save the Robj for the next tool
save(seurat_obj, file = "seurat_obj_subset2.Robj")

# EOF