# TOOL spatial-transcriptomics-seurat-identify-spatial-defined-tissue-domains-HD-v5.R: "Seurat v5 HD -Identify spatially-defined tissue domains" (This tool accurately finds tissue domains rather than cell types of spatial data)
# INPUT seurat_obj_clustering.Robj: "Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL spatiaaliplotti_Spatial.008um.pdf
# OUTPUT OPTIONAL spatiaaliplotti_Spatial.016um.pdf
# OUTPUT OPTIONAL seurat_obj_banksy_Spatial.008um.Robj
# OUTPUT OPTIONAL seurat_obj_banksy_Spatial.016um.Robj
# PARAMETER assay: "Assay to use" TYPE [Spatial.008um: Spatial.008um, Spatial.016um: Spatial.016um] DEFAULT Spatial.008um (Choose between 8 and 16 um bin assays. 8um bin is recommended for analysis)
# PARAMETER OPTIONAL resolution: "Resolution" TYPE DECIMAL DEFAULT 0.5
# PARAMETER OPTIONAL dims.reduction: "dimensions to reduce" TYPE INTEGER DEFAULT 30
# PARAMETER OPTIONAL lazy: "Lazy calculation" TYPE [FALSE, TRUE] DEFAULT FALSE (Makes the analysis faster but no "BANKSY" assay will be created to the seurat object)
# PARAMETER OPTIONAL lambda: "lambda" TYPE DECIMAL DEFAULT 0.5 (A parameter to weight the contributions of the cell-transcriptome matrix and the neighbor expression matrices. Smaller lambda emphasizes cell's own transcriptomes and causes cells to cluster according to cell type. Bigger lambda causes cells to cluster according to tissue domain.)
# PARAMETER OPTIONAL k_geom: "Amount of neighbours" TYPE INTEGER DEFAULT 15 (Local neighborhood size. Larger values will yield larger domains)
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

gc()


# This might become a problem later

# fixes as.graph error
# for (g in Graphs(seurat_obj)) {
# seurat_obj[[g]] <- NULL
# }



# We ended up deciding to use only 1 assay per run, because if lazy = FALSE, banksy creates a new assay called BANKSY based on the assay you choose.
# Thus looping over assays would only take the last assay used likely? Also this is computationally the most expensive step.

DefaultAssay(seurat_obj) <- assay

img_name <- Images(seurat_obj, assay = assay)


# print("Does data slot exist:")
# length(seurat_obj@assays$Spatial.008um@layers$data) > 0

# Maybe an image param has to be created also. Waiting for advice from Meilahti tho

# Should return TRUE. That slot is used in RunBanksy function


if (lazy) {

seurat_obj <- RunBanksy(seurat_obj, lambda = lambda, verbose = TRUE,
    assay = assay, slot = "data",  features = "variable",
    k_geom = k_geom, lazy = T, split.scale = T, parallel = parallel, num_cores = num_cores)

# If lazy = TRUE, no BANKSY assay is created but a reduction called banksy is, just find new neighbors and clusters and plot that


seurat_obj <- FindNeighbors(seurat_obj, reduction = "BANKSY", dims = 1:dims.reduction)
seurat_obj <- FindClusters(seurat_obj, cluster.name = "banksy_cluster", resolution = resolution)


coords <- GetTissueCoordinates(seurat_obj, image = img_name)

common_cells <- intersect(Cells(seurat_obj), rownames(coords))

# for (g in Graphs(seurat_obj)) {
#   seurat_obj[[g]] <- NULL
# }


seurat_obj_sub <- subset(seurat_obj, cells = common_cells)


pdf(file = paste0("spatiaaliplotti_", assay, ".pdf"), width = width, height = height)

    p <- SpatialDimPlot(seurat_obj_sub, images = img_name, group.by = "banksy_cluster", label = T, repel = T, label.size = label.size)
    print(p)


    banksy_cells <- CellsByIdentities(seurat_obj_sub)

    p1 <- SpatialDimPlot(seurat_obj_sub, cells.highlight = banksy_cells[setdiff(names(banksy_cells), "NA")], cols.highlight = c("#FFFF00", "grey50"),
    facet.highlight = T, combine = T)

    print(p1)

dev.off()

seurat_obj <- seurat_obj_sub
rm(seurat_obj_sub)
save(seurat_obj, file = paste0("seurat_obj_banksy_", assay, ".Robj"))


} else {


# Else for lazy = FALSE
# Run Banksy
# The RunBanksy function creates a new BANKSY assay, which can be used for dimensional reduction and clustering:



seurat_obj <- RunBanksy(seurat_obj, lambda = lambda, slot = "data", k_geom = k_geom, lazy = FALSE, verbose = TRUE, assay = DefaultAssay(seurat_obj), parallel = TRUE, num_cores = num_cores)

# Make BANKSY Default assay and run pca, neighboring and clustering
DefaultAssay(seurat_obj) <- "BANKSY"

seurat_obj <- RunPCA(seurat_obj, assay = "BANKSY", reduction.name = "pca.banksy", features = rownames(seurat_obj), npcs = 30)
seurat_obj <- FindNeighbors(seurat_obj, reduction = "pca.banksy", dims = 1:dims.reduction)
seurat_obj <- FindClusters(seurat_obj, resolution = resolution, cluster.name = "banksy_cluster")

# Set banksy_cluster as idents
Idents(seurat_obj) <- "banksy_cluster"

# Plotting

pdf(file = paste0("spatiaaliplotti_", assay, ".pdf"), width = width, height = height)

    p <- SpatialDimPlot(seurat_obj, images = img_name, group.by = "banksy_cluster", label = T, repel = T, label.size = label.size)
    print(p)


    banksy_cells <- CellsByIdentities(seurat_obj)

    p1 <- SpatialDimPlot(seurat_obj, cells.highlight = banksy_cells[setdiff(names(banksy_cells), "NA")], cols.highlight = c("#FFFF00", "grey50"),
    facet.highlight = T, combine = T)

    print(p1)

dev.off()

save(seurat_obj, file = paste0("seurat_obj_banksy_", assay, ".Robj"))

}

# EOF