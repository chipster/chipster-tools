# TOOL spatial-transcriptomics-seurat-clustering-v5-HD.R: "Seurat v5 HD -Clustering" (This tool performs clustering for a single sample or multiple samples that have been combined into one Seurat object. For the 8 µm assay, sketch-based clustering can be enabled: a representative subset of bins is clustered in memory, and the results are then projected back to all bins. This is recommended for large datasets where standard clustering is slow or memory-intensive.)
# INPUT seurat_object.Robj: "Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL seurat_obj_clustering.Robj
# OUTPUT OPTIONAL clustering_plots.pdf
# PARAMETER OPTIONAL dims.reduction: "Dimensions of reduction to use" TYPE INTEGER DEFAULT 30 (Dimensions of reduction to use for clustering and UMAP.)
# PARAMETER OPTIONAL res: "Resolution for granularity for clustering" TYPE DECIMAL DEFAULT 0.8 (Resolution parameter that sets the granularity of the clustering. Higher values lead to greater number of clusters.)
# PARAMETER OPTIONAL use.sketch: "Use sketch-based clustering for 8 µm assay" TYPE [yes: yes, no: no] DEFAULT no (Use sketch-based clustering for the 8 µm assay. A representative subset of bins is selected using leverage scores, clustered in memory, and the cluster labels and UMAP are then projected back to all bins. Recommended for large datasets. The 16 µm assay is always clustered using the standard method.)
# PARAMETER OPTIONAL sketch.ncells: "Number of bins to use in sketch" TYPE INTEGER DEFAULT 50000 (Number of bins to include in the sketch for the 8 µm assay. Only used when sketch-based clustering is enabled. The default of 50000 follows the Seurat v5 Visium HD vignette recommendation.)
# RUNTIME R-4.5.1-seurat5
# SLOTS 4
# TOOLS_BIN ""


# 2026-02 ML
# 2026-06 ML added sketch-based clustering option for 8um assay

# Load seurat object (called seurat_obj)
load("seurat_object.Robj")
library(Seurat)
library(dplyr)
library(Matrix)
library(gplots)
library(ggplot2)
library(patchwork)
library(Biobase)
source(file.path(chipster.common.lib.path, "tool-utils.R"))
#print(package.version("Seurat"))
#documentVersion("Seurat", package.version("Seurat"))
 
assay_names <- Assays(seurat_obj) # e.g. "Spatial.008um" "Spatial.016um"
 
pdf(file = "clustering_plots.pdf", width = 9, height = 12)
 
for (i in 1:length(assay_names)) {
    DefaultAssay(seurat_obj) <- assay_names[i]
    assay_bin   <- assay_names[i]
    just.bin    <- sub("Spatial\\.", "", assay_names[i])
    pca_bin     <- paste("pca.", just.bin, sep = "")
    cluster_bin <- paste("seurat_cluster.", just.bin, sep = "")
    umap_bin    <- paste("umap.", just.bin, sep = "")
 
    # Decide whether to use sketch-based clustering for this assay.
    # Sketching is only applied to the 8um assay when the user has enabled it.
    do.sketch <- (use.sketch == "yes") && grepl("008um", assay_bin)
 
    if (do.sketch) {
        # ---- Sketch-based clustering (8um assay only) ----------------------
        # Follows the Seurat v5 Visium HD vignette exactly.
        # The full dataset stays in memory throughout (no BPCells package needed).
        # Reduction and metadata names match the vignette so that downstream
        # tools can rely on them: cluster labels in "seurat_cluster.projected",
        # full-dataset UMAP in "full.umap.sketch", full PCA in "full.pca.sketch".
 
        # Step 1: build the sketch (leverage-score subsampling).
        # DefaultAssay must be the full 8um assay here.
        # ScaleData is required before SketchData: LeverageScore aborts with
        # "too slow" error if the scaled matrix is absent.
        DefaultAssay(seurat_obj) <- assay_bin
        seurat_obj <- FindVariableFeatures(seurat_obj, verbose = FALSE)
        seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
        seurat_obj <- SketchData(
            object         = seurat_obj,
            ncells         = sketch.ncells,
            method         = "LeverageScore",
            sketched.assay = "sketch",
            features = VariableFeatures(seurat_obj)
        )
 
        # Step 2: cluster the sketch assay
        DefaultAssay(seurat_obj) <- "sketch"
        seurat_obj <- FindVariableFeatures(seurat_obj, verbose = FALSE)
        seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
        seurat_obj <- RunPCA(seurat_obj, assay = "sketch", reduction.name = "pca.sketch", verbose = FALSE)
        seurat_obj <- FindNeighbors(seurat_obj, assay = "sketch", reduction = "pca.sketch", dims = 1:dims.reduction, verbose = FALSE)
        seurat_obj <- FindClusters(seurat_obj, resolution = res, cluster.name = "seurat_cluster.sketched", verbose = FALSE)
        # return.model = TRUE is required so ProjectData can extend the UMAP to all bins
        seurat_obj <- RunUMAP(seurat_obj, reduction = "pca.sketch", reduction.name = "umap.sketch", return.model = TRUE, dims = 1:dims.reduction, verbose = FALSE)
 
        # Step 3: project cluster labels and UMAP back to all 8um bins.
        # After this call:
        #   cluster labels  -> seurat_obj$seurat_cluster.projected  (all bins)
        #   full PCA        -> seurat_obj[["full.pca.sketch"]]
        #   full UMAP       -> seurat_obj[["full.umap.sketch"]]
        seurat_obj <- ProjectData(
            object             = seurat_obj,
            assay              = assay_bin,
            full.reduction     = "full.pca.sketch",
            sketched.assay     = "sketch",
            sketched.reduction = "pca.sketch",
            umap.model         = "umap.sketch",
            dims               = 1:dims.reduction,
            refdata            = list(seurat_cluster.projected = "seurat_cluster.sketched")
        )
        # Copy projected labels into the standard per-assay cluster column so
        # the rest of the loop (SpatialDimPlot) works the same as the non-sketch path.
        seurat_obj[[cluster_bin]] <- seurat_obj[["seurat_cluster.projected"]]
 
        # Switch to the full assay for plotting
        DefaultAssay(seurat_obj) <- assay_bin
 
        # Visualisation: UMAP of all bins (projected) + spatial map
        dim.plot <- DimPlot(seurat_obj, reduction = "full.umap.sketch", group.by = cluster_bin,
                            label = TRUE, repel = TRUE, alpha = 0.1) +
            NoLegend() +
            ggtitle(paste0(just.bin, " (sketch-based, n=", sketch.ncells, ")"))
 
    } else {
        print("Finding neighbors")
        # ---- Standard clustering -------------------------------------------
        seurat_obj <- FindNeighbors(seurat_obj, reduction = pca_bin, dims = 1:dims.reduction, verbose = FALSE)
                print("Finding clusters")

        seurat_obj <- FindClusters(seurat_obj, resolution = res, cluster.name = cluster_bin, verbose = FALSE)
                print("Finding UMAP")

        seurat_obj <- RunUMAP(seurat_obj, reduction = pca_bin, reduction.name = umap_bin, dims = 1:dims.reduction, verbose = FALSE)
 
        dim.plot <- DimPlot(seurat_obj, reduction = umap_bin, group.by = cluster_bin,
                            label = TRUE, repel = TRUE) +
            NoLegend() +
            ggtitle(just.bin)
    }
 
    cluster.plot <- SpatialDimPlot(seurat_obj, group.by = cluster_bin,
                                   pt.size.factor = 1.2, label = TRUE,
                                   repel = TRUE, label.size = 4) +
        theme(legend.position = "right")

    Idents(seurat_obj) <- "seurat_clusters"
    cells <- CellsByIdentities(seurat_obj, idents=c(0,4,30,34,35))
    p <- SpatialDimPlot(seurat_obj, cells.highlight = cells[setdiff(names(cells), "NA")], 
                    cols.highlight = c("#FFFF00","grey50"), facet.highlight = T, combine=T) + NoLegend()
 
    print(dim.plot)
    print(cluster.plot)
    print(p)

} # end loop over assays
 
dev.off()
 
# Save the Robj for the next tool
save(seurat_obj, file = "seurat_obj_clustering.Robj")
 
# EOF
 