# TOOL spatial-transcriptomics-seurat-subset-regions-HD-v5.R: "Seurat v5 HD -Subset regions " (This tool subsets regions based on clusters and or coordinates)
# INPUT seurat_obj_clustering.Robj: "Seurat object" TYPE GENERIC
# INPUT coords.file.csv: "Seurat scRNA data" TYPE GENERIC
# OUTPUT OPTIONAL spatiaaliplotti.pdf
# OUTPUT OPTIONAL spatiaaliplotti2.pdf
# OUTPUT OPTIONAL spatiaaliplotti3.pdf
# PARAMETER coords.file: "Coordination csv file" TYPE [TRUE: yes, FALSE: no] DEFAULT TRUE ()
# PARAMETER coords: "coords..." TYPE [TRUE: yes, FALSE: no] DEFAULT TRUE ()
# PARAMETER OPTIONAL x_coord_min: "Subset x min coordinate" TYPE INTEGER DEFAULT 0
# PARAMETER OPTIONAL x_coord_max: "Subset x max coordinate" TYPE INTEGER DEFAULT 100
# PARAMETER OPTIONAL y_coord_min: "Subset y min coordinate" TYPE INTEGER DEFAULT 0
# PARAMETER OPTIONAL y_coord_max: "Subset y max coordinate" TYPE INTEGER DEFAULT 100
# PARAMETER OPTIONAL label.size: "determine the label size of the plots" TYPE INTEGER DEFAULT 3
# PARAMETER OPTIONAL chosen_clusters: "Subset of clusters" TYPE STRING DEFAULT "1,2,3,4,5" (Clusters to subset. If you list multiple clusters, use comma \(,\) as separator, for example "1,2,3,4".)
# RUNTIME R-4.5.1-seurat5
# SLOTS 6
# TOOLS_BIN ""

# 2026-07-17 JV

chosen_clusters <- as.character(chosen_clusters) # A parameter

chosen_clusters <- trimws(unlist(strsplit(chosen_clusters, ",")))
chosen_clusters <- strtoi(chosen_clusters, base = 0L)

library(Seurat)
library(SeuratObject)
library(ggplot2)
library("sf")

load("seurat_obj_clustering.Robj")



# Make graphs empty because it is causing issues
# seurat_obj <- UpdateSeuratObject(seurat_obj)

#seurat_obj@graphs <- list()

# This loop works, check what causes this graph issue, a possible contact in Meilahti
for (g in Graphs(seurat_obj)) {
seurat_obj[[g]] <- NULL
}


# Recompute graphs if needed
#subset_obj <- RunPCA(seurat_obj, dims = 1:10)
#subset_obj <- FindNeighbors(subset_obj, dims = 1:10)
#subset_obj <- FindClusters(subset_obj)

# Set some clusters as Idents
seurat_obj <- SetIdent(seurat_obj, value = seurat_obj$seurat_clusters)


# Subset 
subset_obj <- subset(seurat_obj, idents = chosen_clusters)

# Check amount of images -> Message sent to Meilahti
print("Images of sobj")
print(names(subset_obj@images))

# The simplest way to do this
pdf(file = "spatiaaliplotti.pdf")
p <- SpatialDimPlot(subset_obj, label = T, crop = T, label.size = 3) + NoLegend() + 
  SpatialDimPlot(subset_obj, label = T, crop = F, label.size = 3)

print(p)

dev.off()


# Maybe add an option to use coords instead or also, then:
if (coords) {

#If coordinates want to be used, then:

# Adds x and y axis numbers for coordinate search
theme_axis_labels <- theme(axis.text.x = element_text(size = 10), 
                           axis.text.y = element_text(size = 10),
                           axis.title.x = element_text(size = 10),
                           axis.title.y = element_text(size = 10))

coord_plot <- SpatialDimPlot(subset_obj, label = T, label.size = 3, crop = T) + theme_axis_labels

x_coord_min # a param that has range: start-end for x axis
x_coord_max # a param that has range: start-end for y axis

y_coord_min
y_coord_max

# Original
# img_name <- Images(subset_obj)[1]
# img_name

# subset_obj[[img_name]] <- Crop(subset_obj[[img_name]], x = c(x_coord_min, x_coord_max), y = c(y_coord_min, y_coord_max))
# }

# For loop for image cropping, may not work
for (img_name in Images(subset_obj)) {
  subset_obj[[img_name]] <- Crop(
    subset_obj[[img_name]],
    x = c(x_coord_min, x_coord_max),
    y = c(y_coord_min, y_coord_max)
  )
}

pdf(file = "spatiaaliplotti2.pdf")

p1 <- SpatialDimPlot(subset_obj, label = T) + theme_axis_labels

print(p1)

dev.off()

}

if (coords.file) {

  img_name <- Images(subset_obj)[1]

  coordinates <- as.data.frame(read.csv("coords.file.csv"))

  head(coordinates)

  segment <- CreateSegmentation(coordinates)

  # Overlay function wants an "sf" package
  seurat_obj[[img_name]] <- Overlay(seurat_obj[[img_name]], segment)
  segment <- subset(seurat_obj, cells = Cells(seurat_obj[[img_name]]))

  pdf("spatiaaliplotti3.pdf")

  p<- SpatialDimPlot(segment)

  print(p)

  dev.off()
}

# EOF
