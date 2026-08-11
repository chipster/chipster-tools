# TOOL spatial-transcriptomics-seurat-identify-spatial-defined-tissue-domains-HD-v5.R: "Seurat v5 HD -Identify spatially-defined tissue domains" (This tool identifies differentially expressed genes between two user defined clusters and visualizes these genes on top of the tissue image.)
# INPUT seurat_obj_clustering.Robj: "Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL spatiaaliplotti.pdf
# OUTPUT OPTIONAL seurat_obj_banksy.Robj
# PARAMETER OPTIONAL resolution: "reso" TYPE DECIMAL DEFAULT 0.5
# PARAMETER OPTIONAL dims.reduction: "dims to redu" TYPE INTEGER DEFAULT 30
# PARAMETER OPTIONAL lazy: "Lazy calc or not" TYPE [FALSE, TRUE] DEFAULT FALSE
# PARAMETER OPTIONAL lambda: "lambda for weighting" TYPE DECIMAL DEFAULT 0.5
# PARAMETER OPTIONAL k_geom: "Amount of neighbours" TYPE INTEGER DEFAULT 15
# PARAMETER OPTIONAL label.size: "determine the label size of the plots" TYPE INTEGER DEFAULT 3
# PARAMETER OPTIONAL width: "Width of the pdf" TYPE INTEGER DEFAULT 10
# PARAMETER OPTIONAL height: "Height of the pdf" TYPE INTEGER DEFAULT 10
# RUNTIME R-4.5.1-seurat5
# SLOTS 10
# TOOLS_BIN ""


resolution <- as.numeric(resolution)
dims.reduction <- as.numeric(dims.reduction)
lambda <- as.numeric(lambda)
k_geom <- as.numeric(k_geom)
label.size <- as.numeric(label.size)
lazy <- as.logical(lazy)



num_cores = as.numeric(chipster.threads.max)
parallel = TRUE

# 2026-07-15 JV

# Matrix version 1.8.0 needed for zgeMatrix. Trying to install it on the VM Runtime R-4.5.1-seurat5 through github
library("Seurat")
library("Banksy")
library("SeuratWrappers")

# Current setup: Banksy 1.9.2, matrix 1.7.5, seuratwrappers 0.4.0

# sessionInfo() shows all versions if needed
set.seed(123)

load("seurat_obj_clustering.Robj")

# For testing, make seurat smaller:

# n_cells <- 50000  # target number of cells
# set.seed(42)      # for reproducibility

# cells_to_keep <- sample(colnames(seurat_obj), size = n_cells)
# seurat_obj <- subset(seurat_obj, cells = cells_to_keep)

gc()
# fixes as.graph error

for (g in Graphs(seurat_obj)) {
seurat_obj[[g]] <- NULL
}

# Spatial.008um
# Maybe add this as a param also

print("Assays in seurat obj")
Assays(seurat_obj)

# Loop over the first two assays (8 and 16um bins, not sketch!). If it reaches sketch, dim error...

for (assay in Assays(seurat_obj)[1:2]) {

DefaultAssay(seurat_obj) <- assay

img_name <- Images(seurat_obj, assay = assay)
print("image name")
print(img_name)

# print("Does data slot exist:")
# length(seurat_obj@assays$Spatial.008um@layers$data) > 0

# Maybe an image param has to be created also. Waiting for advice from Meilahti tho

# Should return TRUE. That slot is used in RunBanksy function


if (lazy) {

seurat_obj <- RunBanksy(seurat_obj, lambda = lambda, verbose = TRUE,
    assay = assay, slot = "data",  features = "variable",
    k_geom = k_geom, lazy = T, split.scale = T, parallel = parallel, num_cores = num_cores)

# If lazy = TRUE, no BANKSY assay is created but a reduction called banksy is, just find new neighbors and clusters and plot that
# Like 10000000x faster

print("Banksy done, gc")

seurat_obj <- FindNeighbors(seurat_obj, reduction = "BANKSY", dims = 1:dims.reduction)
print("Neighbor done")
seurat_obj <- FindClusters(seurat_obj, cluster.name = "banksy_cluster", resolution = resolution)
print("Clusters done")


print("Starting coords")
coords <- GetTissueCoordinates(seurat_obj, image = img_name)
print("Coords done")
common_cells <- intersect(Cells(seurat_obj), rownames(coords))

print("Common cells done")
for (g in Graphs(seurat_obj)) {
  seurat_obj[[g]] <- NULL
}

print("For loop done")

seurat_obj_sub <- subset(seurat_obj, cells = common_cells)


print("Subset done, next plotting:")
pdf(file = paste0("spatiaaliplotti_", assay, ".pdf"), width = width, height = height)

    p <- SpatialDimPlot(seurat_obj_sub, images = img_name, group.by = "banksy_cluster", label = T, repel = T, label.size = label.size)
    print(p)


    banksy_cells <- CellsByIdentities(seurat_obj_sub)

    p1 <- SpatialDimPlot(seurat_obj_sub, cells.highlight = banksy_cells[setdiff(names(banksy_cells), "NA")], cols.highlight = c("#FFFF00", "grey50"),
    facet.highlight = T, combine = T)

    print(p1)

dev.off()

save(seurat_obj_sub, file = paste0("seurat_obj_banksy_", assay, ".Robj"))

print("Loop finished, new round!!!")

}
}

# } else {


# # Else for lazy = FALSE, commented out for test purposes.
# # Run Banksy
# seurat_obj <- RunBanksy(seurat_obj, lambda = lambda, slot = "data", k_geom = k_geom, lazy = FALSE, verbose = TRUE, assay = DefaultAssay(seurat_obj), parallel = TRUE, num_cores = num_cores)

# # Make BANKSY Default assay and run pca, neighboring and clustering
# DefaultAssay(seurat_obj) <- "BANKSY"

# print("PCA: ")
# seurat_obj <- RunPCA(seurat_obj, assay = "BANKSY", reduction.name = "pca.banksy", features = rownames(seurat_obj), npcs = 30)

# print("Neighbours")
# seurat_obj <- FindNeighbors(seurat_obj, reduction = "pca.banksy", dims = 1:dims.reduction)

# print("Clusters: ")
# seurat_obj <- FindClusters(seurat_obj, resolution = resolution, cluster.name = "banksy_cluster")

# # Set banksy_cluster as idents
# Idents(seurat_obj) <- "banksy_cluster"

# # Plotting

# pdf(file = paste0("spatiaaliplotti_", assay, ".pdf"), width = width, height = height)

#     p <- SpatialDimPlot(seurat_obj, images = img_name, group.by = "banksy_cluster", label = T, repel = T, label.size = label.size)
#     print(p)


#     banksy_cells <- CellsByIdentities(seurat_obj)

#     p1 <- SpatialDimPlot(seurat_obj, cells.highlight = banksy_cells[setdiff(names(banksy_cells), "NA")], cols.highlight = c("#FFFF00", "grey50"),
#     facet.highlight = T, combine = T)

#     print(p1)

# dev.off()

# save(seurat_obj, file = paste0("seurat_obj_banksy_", assay, ".Robj"))

# }

# EOF
