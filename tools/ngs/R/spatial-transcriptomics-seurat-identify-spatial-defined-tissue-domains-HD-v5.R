# TOOL spatial-transcriptomics-seurat-identify-spatial-defined-tissue-domains-HD-v5.R: "Seurat v5 -Identify spatially-defined tissue domains with BANKSY" (This tool uses BANKSY to incorporate neighborhood information for clustering of Visium HD data. BANKSY augments the expression pattern of a spot with both the mean and the gradient of gene expression levels in a broader neighborhood of a spot.)
# INPUT seurat_obj_clustering.Robj: "Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL seurat_obj_banksy_Spatial.008um.Robj
# OUTPUT OPTIONAL seurat_obj_banksy_Spatial.016um.Robj
# OUTPUT OPTIONAL BANKSY_plots.pdf
# PARAMETER assay: "Assay to use" TYPE [Spatial.008um: Spatial.008um, Spatial.016um: Spatial.016um] DEFAULT Spatial.008um (Choose between 8 and 16 um bin assays. 8um bin is recommended for analysis)
# PARAMETER OPTIONAL lambda: "Lambda" TYPE DECIMAL FROM 0 TO 1 DEFAULT 0.8 (Influence of the neighbourhood. Bigger lambda causes spots to cluster according to tissue domain, while smaller lambda emphasizes spot's own transcriptomes and causes spots to cluster according to cell type.)
# PARAMETER OPTIONAL k_geom: "Number of neighbors" TYPE INTEGER DEFAULT 50 (Local neighborhood size. This number of spots defines a spot’s neighborhood. Larger values will yield larger domains.)
# PARAMETER OPTIONAL resolution: "Resolution value for clustering" TYPE DECIMAL DEFAULT 0.8
# PARAMETER OPTIONAL dims.reduction: "Number of principal components to use for clustering" TYPE INTEGER DEFAULT 30
# PARAMETER OPTIONAL label.size: "Label size of the plots" TYPE INTEGER DEFAULT 3
# PARAMETER OPTIONAL width: "Width of the plot file" TYPE INTEGER DEFAULT 10
# PARAMETER OPTIONAL height: "Height of the plot file" TYPE INTEGER DEFAULT 10
# RUNTIME R-4.5.1-visium-hd
# SLOTS 10
# TOOLS_BIN ""

# For png, cairo should be installed in this RUNTIME
# DONE -- BANKSY also working now

resolution <- as.numeric(resolution)
dims.reduction <- as.numeric(dims.reduction)
lambda <- as.numeric(lambda)
k_geom <- as.numeric(k_geom)
label.size <- as.numeric(label.size)

# For now, use lazy = TRUE as default

## PARAMETER OPTIONAL lazy: "Lazy calculation" TYPE [FALSE, TRUE] DEFAULT TRUE (Makes the analysis faster, but no "BANKSY" assay will be created to the Seurat object.)

lazy <- TRUE

# For now, use only variable features for computing BANKSY clusters. Banksy does allow all genes to be used, but this is computationally expensive and not needed for most cases.
# PARAMETER OPTIONAL features: "Which genes to use for computing tissue domains" TYPE [all: "all", variable: "variable only"] DEFAULT variable (Select either variable or all genes for computing tissue domains. Only applies if Lazy calculation is enabled.)
features <- as.character("variable")
num_cores <- as.numeric(chipster.threads.max)
parallel <- TRUE

# 2026-07-15 JV
# 2026-10-15 JV Make lazy = TRUE always & use only variable features
# 2026-10-15 JV Add png solution to make pdf smaller in size

library("Seurat")
library("Banksy")
library("SeuratWrappers")
library("ragg")
library("magick")

# Current setup: Banksy 1.9.2, matrix 1.7.5, seuratwrappers 0.4.0

load("seurat_obj_clustering.Robj")

# This might become a problem later
# fixes as.graph error
# for (g in Graphs(seurat_obj)) {
# seurat_obj[[g]] <- NULL
# }

# We ended up deciding to use only 1 assay per run, because if lazy = FALSE, banksy creates a new assay called BANKSY based on the assay you have chosen.
# Thus looping over assays would only take the last assay used likely? Also, this is computationally the most expensive step.

DefaultAssay(seurat_obj) <- assay


# Image is needed for SpatialDimPlot, otherwise an error pops out (specify the image based on your assay)
img_name <- Images(seurat_obj, assay = assay)


if (lazy) {
    seurat_obj <- RunBanksy(
        seurat_obj,
        lambda = lambda,
        verbose = TRUE,
        assay = assay,
        slot = "data",
        features = features,
        k_geom = k_geom,
        lazy = TRUE,
        split.scale = TRUE,
        parallel = parallel,
        num_cores = num_cores
    )

    # If lazy = TRUE, no BANKSY assay is created, but a reduction called banksy is, so we find new neighbors and clusters and plot that.

    seurat_obj <- FindNeighbors(seurat_obj, reduction = "BANKSY", dims = 1:dims.reduction)
    seurat_obj <- FindClusters(seurat_obj, cluster.name = "banksy_cluster", resolution = resolution)

    Idents(seurat_obj) <- "banksy_cluster"

    #   pdf(file = paste0("BANKSY_plot_", assay, ".pdf"), width = width, height = height)

    #   p <- SpatialDimPlot(seurat_obj, images = img_name, group.by = "banksy_cluster", label = T, repel = T, label.size = label.size)
    #   print(p)

    #   banksy_cells <- CellsByIdentities(seurat_obj)

    #   p1 <- SpatialDimPlot(seurat_obj, images = img_name, cells.highlight = banksy_cells[setdiff(names(banksy_cells), "NA")], cols.highlight = c("#FFFF00", "grey50"),
    #                        facet.highlight = T, combine = T)

    #   print(p1)

    # Leave them out for now, these are the huge > 500Mb pdfs

    # OUTPUT OPTIONAL BANKSY_plot_Spatial.008um.pdf
    # OUTPUT OPTIONAL BANKSY_plot_Spatial.016um.pdf

    #   dev.off()

    # Output directory for PNG plots

    output_dir <- paste0(getwd(), "/Banksy_plots")
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

    for (i in seq_along(assay)) {
        assay_bin <- assay[i]
        DefaultAssay(seurat_obj) <- assay_bin
        just.bin <- sub("Spatial\\.", "", assay_bin)

        p0 <- SpatialDimPlot(
            seurat_obj,
            images = img_name,
            group.by = "banksy_cluster",
            label = T,
            repel = T,
            label.size = label.size
        )

        banksy_cells <- CellsByIdentities(seurat_obj)

        p1 <- SpatialDimPlot(
            seurat_obj,
            images = img_name,
            cells.highlight = banksy_cells[setdiff(names(banksy_cells), "NA")],
            cols.highlight = c("#FFFF00", "grey50"),
            facet.highlight = T,
            combine = T
        )

        # Put plots into a list
        plots <- list(p0, p1)

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

    image_write(images, format = "pdf", path = "BANKSY_plots.pdf")

    save(seurat_obj, file = paste0("seurat_obj_banksy_", assay, ".Robj"))
}

# else {

#   # Else for lazy = FALSE
#   # The RunBanksy function creates a new BANKSY assay, which can be used for dimensional reduction and clustering:
#   # features parameter only works for lazy=T because RunBanksy() has hardcoded limitations to sizes. "all" features would be too much, at least for this dummy data.

#   seurat_obj <- RunBanksy(seurat_obj, lambda = lambda, slot = "data", features = "variable", k_geom = k_geom, lazy = FALSE, verbose = TRUE, assay = assay, parallel = TRUE, num_cores = num_cores)

#   # Make BANKSY Default assay and run pca, neighboring and clustering
#   DefaultAssay(seurat_obj) <- "BANKSY"

#   seurat_obj <- RunPCA(seurat_obj, assay = "BANKSY", reduction.name = "pca.banksy", features = rownames(seurat_obj), npcs = 30)
#   seurat_obj <- FindNeighbors(seurat_obj, reduction = "pca.banksy", dims = 1:dims.reduction)
#   seurat_obj <- FindClusters(seurat_obj, resolution = resolution, cluster.name = "banksy_cluster")

#   # Set banksy_cluster as idents
#   Idents(seurat_obj) <- "banksy_cluster"

#   # Plotting

#   # pdf(file = paste0("BANKSY_plot_", assay, ".pdf"), width = width, height = height)

#   # p <- SpatialDimPlot(seurat_obj, images = img_name, group.by = "banksy_cluster", label = T, repel = T, label.size = label.size)
#   # print(p)

#   # banksy_cells <- CellsByIdentities(seurat_obj)

#   # p1 <- SpatialDimPlot(seurat_obj, images = img_name, cells.highlight = banksy_cells[setdiff(names(banksy_cells), "NA")], cols.highlight = c("#FFFF00", "grey50"),
#   #                      facet.highlight = T, combine = T)

#   # print(p1)

#   # dev.off()

#   banksy_cells <- CellsByIdentities(seurat_obj)

# # Output directory for PNG plots

# output_dir <- paste0(getwd(), "/Banksy_plots")
# dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# for (i in seq_along(assay)) {

#   assay_bin <- assay[i]
#   DefaultAssay(seurat_obj) <- assay_bin
#   just.bin <- sub("Spatial\\.", "", assay_bin)

#   p0 <- SpatialDimPlot(seurat_obj, images = img_name, group.by = "banksy_cluster", label = T, repel = T, label.size = label.size)

#   banksy_cells <- CellsByIdentities(seurat_obj)

#   p1 <- SpatialDimPlot(seurat_obj, images = img_name, cells.highlight = banksy_cells[setdiff(names(banksy_cells), "NA")], cols.highlight = c("#FFFF00", "grey50"),
#                        facet.highlight = T, combine = T)

#   # Put plots into a list
#   plots <- list(p0, p1)

#   # Loop over the plots and save each as a separate PNG
#   for (j in seq_along(plots)) {

#     filename <- file.path(
#       output_dir,
#       paste0(just.bin, "_plot_", j, ".png")
#     )

#     agg_png(
#       filename = filename,
#       width = 10,
#       height = 10,
#       res = 300,
#       bitsize = 16,
#       units = "in"
#     )

#     print(plots[[j]])

#     dev.off()
#   }
# }

# png_files <- list.files(
#   output_dir,
#   pattern = "\\.png$",
#   full.names = TRUE
# )

# # Sort them so the PDF follows assay -> plot order
# png_files <- sort(png_files)

# # Read all PNGs
# images <- image_read(png_files)

# image_write(images, format="pdf", path = "BANKSY_plots.pdf")

#   save(seurat_obj, file = paste0("seurat_obj_banksy_", assay, ".Robj"))

# }

# EOF
