# TOOL spatial-transcriptomics-seurat-integration-with-scRNA-HD-v5.R: "Seurat v5 HD -Integration with scRNA-seq data" (This tool identifies differentially expressed genes between two user defined clusters and visualizes these genes on top of the tissue image.)
# INPUT seurat_obj_clustering.Robj: "Seurat object" TYPE GENERIC
# INPUT scRNAseq_ref.Rds: "Seurat scRNA data" TYPE GENERIC
# OUTPUT OPTIONAL spatiaaliplotti.pdf
# PARAMETER assay: "Assay to use for RCTD" TYPE [Spatial.008um: "Spatial.008um", Spatial.016um: "Spatial.016um"] DEFAULT Spatial.008um
# PARAMETER OPTIONAL label.size: "determine the label size of the plots" TYPE INTEGER DEFAULT 3
# RUNTIME R-4.5.1-seurat5
# SLOTS 20
# TOOLS_BIN ""

assay <- as.character(assay)
resolution = 0.5
dims.reduction = 30
label.size = 2
num_cores = as.numeric(chipster.threads.max)
parallel = TRUE

# 2026-07-15 JV


# Matrix version 1.8.0 needed for zgeMatrix. Trying to install it on the VM Runtime R-4.5.1-seurat5 through github

library(Seurat)
library(ggplot2)
library("Matrix")
library("SeuratWrappers")
library("spacexr")


options(future.globals.maxSize = 3000 * 1024^2)  # 3000 MiB
#print("Session info:")
#sessionInfo()

load("seurat_obj_clustering.Robj")
# n_cells <- 20000  # target number of cells
# set.seed(42)      # for reproducibility

# cells_to_keep <- sample(colnames(seurat_obj), size = n_cells)
# seurat_obj <- subset(seurat_obj, cells = cells_to_keep)

spatial_obj <- seurat_obj
rm(seurat_obj)



print("Loading ref object")

#load("seurat_obj_scrna.Robj")


# cells_to_keep <- sample(colnames(data), size = n_cells)
# data <- subset(data, cells = cells_to_keep)

ref <- readRDS("scRNAseq_ref.Rds")
#rm(data)


# Basic filter to make sure
ref <- subset(ref, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)


Idents(ref) <- "subclass_label"
counts <- ref[["RNA"]]$counts
cluster <- as.factor(ref$subclass_label)
nUMI <- ref$nCount_RNA
levels(cluster) <- gsub("/", "-", levels(cluster))
cluster <- droplevels(cluster)

# create the RCTD reference object
reference <- Reference(counts, cluster, nUMI)


if ("sketch" %in% Assays(spatial_obj)) {

print("Sketch in the assay was found, using that:")

counts_hd <- spatial_obj[["sketch"]]$counts
spatial_obj_cells_hd <- colnames(spatial_obj[["sketch"]])

# Claude help
coords <- GetTissueCoordinates(spatial_obj)[spatial_obj_cells_hd,1:2]
rownames(coords) <- coords$cell
coords <- coords[spatial_obj_cells_hd,1:2]

# create the RCTD query object
query <- SpatialRNA(coords, counts_hd, colSums(counts_hd))

# Run RCTD
print("Starting RCCTD")
RCTD <- create.RCTD(query, reference, max_cores = num_cores, CELL_MIN_INSTANCE = 0)
RCTD <- run.RCTD(RCTD)

spatial_obj <- AddMetaData(spatial_obj, metadata = RCTD@results$results_df)

print("RCTDD DONE")
# Project RCTD labels 

# project RCTD labels from sketched cortical cells to all cortical cells
spatial_obj$first_type <- as.character(spatial_obj$first_type)
spatial_obj$first_type[is.na(spatial_obj$first_type)] <- 'Unknown'
spatial_obj <- ProjectData(
  object = spatial_obj,
  assay = "Spatial.008um",
  full.reduction = "pca.spatial_obj",
  sketched.assay = "sketch",
  sketched.reduction = "pca.spatial_obj.sketch",
  umap.model = "umap.spatial_obj.sketch",
  dims = 1:50,
  refdata = list(full_first_type = "first_type")
)

# Plotting


DefaultAssay(spatial_obj) <- "Spatial.008um"
print("Final things")
# we only ran RCTD on the cortical cells
# set labels to all other cells as "Unknown"
spatial_obj$full_first_type[is.na(spatial_obj$full_first_type)] <- "Unknown"
Idents(spatial_obj) <- 'full_first_type'

# now we can spatially map the location of any scRNA-seq cell type
# start with Layered (starts with L), excitatory neurons in the cortex
cells <- CellsByIdentities(spatial_obj)
excitatory_names <- sort(grep("^L.* CTX",names(cells),value = TRUE)) # This needs to be customized


pdf(file = "spatiaaliplotti.pdf")
p <- SpatialDimPlot(spatial_obj, cells.highlight = cells[excitatory_names], cols.highlight = c("#FFFF00","grey50"), facet.highlight = T, combine=T, ncol=4)
print(p)

dev.off()

} else {

print(paste0("Sketch in the assay was NOT found, using: ", assay))
print(assay)


DefaultAssay(spatial_obj) <- assay

counts_hd <- spatial_obj[[assay]]$counts
spatial_obj_cells_hd <- colnames(spatial_obj[[assay]])


coords <- GetTissueCoordinates(spatial_obj)
rownames(coords) <- coords$cell
coords <- coords[spatial_obj_cells_hd, c("x", "y")]

# create the RCTD query object
query <- SpatialRNA(coords, counts_hd, colSums(counts_hd))

# Run RCTD
print("Starting RCCTD")
suppressMessages({
RCTD <- create.RCTD(query, reference, max_cores = num_cores, CELL_MIN_INSTANCE = 0)
RCTD <- run.RCTD(RCTD)

spatial_obj <- AddMetaData(spatial_obj, metadata = RCTD@results$results_df)

print("RCTDD DONE")
# Project RCTD labels 

# project RCTD labels from sketched cortical cells to all cortical cells
spatial_obj$first_type <- as.character(spatial_obj$first_type)
spatial_obj$first_type[is.na(spatial_obj$first_type)] <- 'Unknown'


# Plotting


DefaultAssay(spatial_obj) <- assay
print("Final things")
# we only ran RCTD on the cortical cells
# set labels to all other cells as "Unknown"
spatial_obj$full_first_type[is.na(spatial_obj$full_first_type)] <- "Unknown"
Idents(spatial_obj) <- 'full_first_type'

# now we can spatially map the location of any scRNA-seq cell type
# start with Layered (starts with L), excitatory neurons in the cortex
cells <- CellsByIdentities(spatial_obj)
excitatory_names <- sort(grep("^L.* CTX",names(cells),value = TRUE)) # This needs to be customized


pdf(file = "spatiaaliplotti.pdf")
p <- SpatialDimPlot(spatial_obj, cells.highlight = cells[excitatory_names], cols.highlight = c("#FFFF00","grey50"), facet.highlight = T, combine=T, ncol=4)
print(p)

dev.off()
})

}
# EOF






# # Seurat vignette



# #sketch the cortical subset of the Visium HD dataset
# DefaultAssay(seurat_obj) <- "Spatial.008um"
# seurat_obj <- FindVariableFeatures(seurat_obj)

# if (!("sketch" %in% Assays(seurat_obj))) {

# seurat_obj <- SketchData(
#   object = seurat_obj,
#   ncells = 50000,
#   method = "LeverageScore",
#   sketched.assay = "sketch"
# )

# DefaultAssay(seurat_obj) <- "sketch"
# seurat_obj <- ScaleData(seurat_obj)
# seurat_obj <- RunPCA(seurat_obj, assay="sketch", reduction.name = "pca.seurat_obj.sketch", verbose = T)
# seurat_obj <- FindNeighbors(seurat_obj, reduction = "pca.seurat_obj.sketch", dims = 1:50)
# seurat_obj <- RunUMAP(seurat_obj, reduction = "pca.seurat_obj.sketch", reduction.name = "umap.seurat_obj.sketch", return.model = T, dims = 1:50, verbose = T)

# }


# # load in the reference scRNA-seq dataset
# ref <- readRDS("scRNAseq_ref.Rds")


# Idents(ref) <- "subclass_label"
# counts <- ref[["RNA"]]$counts
# cluster <- as.factor(ref$subclass_label)
# nUMI <- ref$nCount_RNA
# levels(cluster) <- gsub("/", "-", levels(cluster))
# cluster <- droplevels(cluster)

# # create the RCTD reference object
# reference <- Reference(counts, cluster, nUMI)

# counts_hd <- seurat_obj[["sketch"]]$counts
# seurat_obj_cells_hd <- colnames(seurat_obj[["sketch"]])
# coords <- GetTissueCoordinates(seurat_obj)[seurat_obj_cells_hd,1:2]

# # create the RCTD query object
# query <- SpatialRNA(coords, counts_hd, colSums(counts_hd))



# # run RCTD
# RCTD <- create.RCTD(query, reference, max_cores = num_cores)
# RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
# # add results back to Seurat object
# seurat_obj <- AddMetaData(seurat_obj, metadata = RCTD@results$results_df)




# # project RCTD labels from sketched cortical cells to all cortical cells
# seurat_obj$first_type <- as.character(seurat_obj$first_type)
# seurat_obj$first_type[is.na(seurat_obj$first_type)] <- 'Unknown'
# seurat_obj <- ProjectData(
#   object = seurat_obj,
#   assay = "Spatial.008um",
#   full.reduction = "pca.seurat_obj.sketch",
#   sketched.assay = "sketch",
#   sketched.reduction = "pca.seurat_obj.sketch",
#   umap.model = "umap.seurat_obj.sketch",
#   dims = 1:50,
#   refdata = list(full_first_type = "first_type")
# )




# DefaultAssay(seurat_obj) <- "Spatial.008um"

# # we only ran RCTD on the cortical cells
# # set labels to all other cells as "Unknown"
# seurat_obj[[]][, "full_first_type"] <- "Unknown"
# seurat_obj$full_first_type[Cells(seurat_obj)] <- seurat_obj$full_first_type[Cells(seurat_obj)]
# Idents(seurat_obj) <- 'full_first_type'

# # now we can spatially map the location of any scRNA-seq cell type
# # start with Layered (starts with L), excitatory neurons in the cortex
# cells <- CellsByIdentities(seurat_obj)
# excitatory_names <- sort(grep("^L.* CTX",names(cells),value = TRUE))

# pdf(file = "spatiaaliplotti.pdf")

# p <- SpatialDimPlot(seurat_obj, cells.highlight = cells[excitatory_names], cols.highlight = c("#FFFF00","grey50"), facet.highlight = T, combine=T, ncol=4)
# print(p)





# # New


# plot_cell_types <- function(data, label) {
#   p <- ggplot(data, aes(x = get(label), y = n, fill = full_first_type)) +
#     geom_bar(stat = "identity", position = "stack") +
#     geom_text(aes(label = ifelse(n >= min_count_to_show_label, full_first_type, "")), position = position_stack(vjust = 0.5), size = 2) +
#     xlab(label) +
#     ylab("# of Spots") +
#     ggtitle(paste0("Distribution of Cell Types across ", label)) +
#     theme_minimal()
# }

# cell_type_banksy_counts <- seurat_obj[[]] %>%
#   dplyr::filter(full_first_type %in% excitatory_names) %>%
#   dplyr::count(full_first_type, banksy_cluster)

# min_count_to_show_label <- 20

# p1 <- plot_cell_types(cell_type_banksy_counts, "banksy_cluster")
# print(p1)




# Idents(seurat_obj) <- 'banksy_cluster'
# seurat_obj$layer_id <- 'Unknown'
# seurat_obj$layer_id[WhichCells(seurat_obj,idents = c(5))] <- "Layer 2/3"
# seurat_obj$layer_id[WhichCells(seurat_obj,idents = c(12))] <- "Layer 4"
# seurat_obj$layer_id[WhichCells(seurat_obj,idents = c(7))] <- "Layer 5"
# seurat_obj$layer_id[WhichCells(seurat_obj,idents = c(3))] <- "Layer 6"






# # set ID to RCTD label
# Idents(seurat_obj) <- 'full_first_type'

# # Visualize distribution of 4 interneuron subtypes
# inhibitory_names <- c("Sst","Pvalb","Vip","Lamp5")
# cell_ids <- CellsByIdentities(seurat_obj, idents = inhibitory_names)
# p2 <- SpatialDimPlot(seurat_obj, cells.highlight = cell_ids, cols.highlight = c("#FFFF00","grey50"), facet.highlight = T, combine=T, ncol=4)
# print(p2)







# # create barplot to show proportions of cell types of interest
# layer_table <- table(seurat_obj$full_first_type, seurat_obj$layer_id)[inhibitory_names,1:4]

# neuron_props <- reshape2::melt(prop.table(layer_table), margin = 1)
# p3 <- ggplot(neuron_props, aes(x = Var1, y = value, fill = Var2)) +
#   geom_bar(stat = "identity", position = "fill") +
#   labs(x = "Cell type", y = "Proportion", fill = "Layer") +
#   theme_classic()
# print(p3)

# dev.off()

# # EOF