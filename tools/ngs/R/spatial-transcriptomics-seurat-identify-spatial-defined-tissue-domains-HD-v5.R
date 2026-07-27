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
# RUNTIME R-4.5.1-seurat5
# SLOTS 6
# TOOLS_BIN ""


resolution <- as.numeric(resolution)
dims.reduction <- as.numeric(dims.reduction)
lambda <- as.numeric(lambda)
k_geom <- as.numeric(k_geom)
label.size <- as.numeric(label.size)
lazy <- as.logical(lazy)

# Chipster performance
# Multithreading and parallelization

num_cores = as.numeric(chipster.threads.max)
parallel = TRUE

# 2026-07-15 JV

# Matrix version 1.8.0 needed for zgeMatrix. Trying to install it on the VM Runtime R-4.5.1-seurat5 through github
library("Seurat")
library("Banksy")
library("SeuratWrappers")

# Current setup: Banksy 1.9.2, matrix 1.7.5, seuratwrappers 0.4.0

# sessionInfo() shows all versions if needed

load("seurat_obj_clustering.Robj")

# Spatial.008um
# Maybe add this as a param also
DefaultAssay(seurat_obj) <- "Spatial.008um"

# Run Banksy
seurat_obj <- RunBanksy(seurat_obj, lambda = 0.8, k_geom = 5, lazy = lazy, verbose = FALSE, assay = DefaultAssay(seurat_obj), parallel = parallel, num_cores = num_cores)

# Make BANKSY Default assay and run pca, neighboring and clustering
DefaultAssay(seurat_obj) <- "BANKSY"

print("PCA: ")
seurat_obj <- RunPCA(seurat_obj, assay = "BANKSY", reduction.name = "pca.banksy", features = rownames(seurat_obj), npcs = 30)

print("Neighbours")
seurat_obj <- FindNeighbors(seurat_obj, reduction = "pca.banksy", dims = 1:dims.reduction)

print("Clusters: ")
seurat_obj <- FindClusters(seurat_obj, resolution = resolution, cluster.name = "banksy_cluster")

# Set banksy_cluster as idents
Idents(seurat_obj) <- "banksy_cluster"

# Plotting

pdf(file = "spatiaaliplotti.pdf")

    p <- SpatialDimPlot(seurat_obj, images = "slice1.008um", group.by = "banksy_cluster", label = T, repel = T, label.size = label.size)
    print(p)


    banksy_cells <- CellsByIdentities(seurat_obj)

    p1 <- SpatialDimPlot(seurat_obj, cells.highlight = banksy_cells[setdiff(names(banksy_cells), "NA")], cols.highlight = c("#FFFF00", "grey50"),
    facet.highlight = T, combine = T)

    print(p1)

dev.off()

save(seurat_obj, file = "seurat_obj_banksy.Robj")

# EOF
