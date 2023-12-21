# TOOL single-cell-seurat-merge-samples.R: "Seurat v5 -Merge samples" (Merge multiple samples for joined analysis.)
# INPUT samples{...}.Robj: "Samples to combine" TYPE GENERIC
# OUTPUT seurat_obj_merged.Robj
# OUTPUT OPTIONAL Dispersion_plot.pdf
# OUTPUT OPTIONAL PCAplots.pdf
# OUTPUT OPTIONAL PCAloadings.txt
# PARAMETER OPTIONAL normalisation.method: "Normalization method to perform" TYPE [LogNormalize:"Global scaling normalization", SCT:"SCTransform"] DEFAULT LogNormalize (Normalize data with global scaling normalization or SCTransform.)
# PARAMETER OPTIONAL totalexpr: "Scaling factor in the normalization" TYPE INTEGER DEFAULT 10000 (Scale each cell to this total number of transcripts.)
# PARAMETER OPTIONAL num.features: "Number of variable genes to return" TYPE INTEGER DEFAULT 2000 (Number of features to select as top variable features, i.e. how many features returned.)
# PARAMETER OPTIONAL num.of.pcas: "Number of PCs to compute" TYPE INTEGER DEFAULT 30 (How many principal components to compute and store. If you get an error message, try lowering the number. This might happen especially if you have low number of cells in your data.)
# PARAMETER OPTIONAL num.of.heatmaps: "Number of principal components to plot as heatmaps" TYPE INTEGER DEFAULT 12 (How many principal components to plot as heatmaps.)
# PARAMETER OPTIONAL loadings: "Print loadings in a file" TYPE [TRUE: yes, FALSE: no] DEFAULT FALSE (Print the PC loadings to a txt file.)
# PARAMETER OPTIONAL num.of.genes.loadings: "Number of genes to list in the loadings file" TYPE INTEGER DEFAULT 5 (How many genes to list in the loadings txt file.)
# PARAMETER OPTIONAL filter.cell.cycle: "Regress out cell cycle differences" TYPE [no:no, all.diff:"all differences", diff.phases:"the difference between the G2M and S phase scores"] DEFAULT no (Would you like to regress out cell cycle scores during data scaling? If yes, should all signal associated with cell cycle be removed, or only the difference between the G2M and S phase scores.)
# RUNTIME R-4.3.2-single-cell
# SLOTS 3
# TOOLS_BIN ""

# 2023-11-29 IH

library(Seurat)
library(dplyr)
library(Matrix)
library(gplots)
library(ggplot2)
options(Seurat.object.assay.version = "v5")

# Load the R-Seurat-objects (called seurat_obj -that's why we need to rename them here)
input.names <- read.table("chipster-inputs.tsv", header = FALSE, sep = "\t")
for (i in 1:nrow(input.names)) {
    load(input.names[i, 1])
    name.of.obj <- paste("seurat_obj_", i, sep = "")
    assign(name.of.obj, seurat_obj)
}
seurat.objects.list <- as.list(mget(objects(pattern = "seurat_obj_")))

# merge first seurat object in the list with the other seurat object(s)
seurat_obj <- merge(x = seurat.objects.list[[1]], y = seurat.objects.list[-1])
seurat_obj

# For cell cycle filtering, read in a list of cell cycle markers, from Tirosh et al, 2015
# cell cycle file is in the container image, not in tools-bin /opt/chipster/tools vs. /opt/chipster/tools-bin
cc.genes <- readLines(con = file.path("/opt/chipster/tools/seurat/regev_lab_cell_cycle_genes.txt"))
# cc.genes <- readLines(con = file.path(chipster.tools.path, "seurat/regev_lab_cell_cycle_genes.txt"))

# We can segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]

# Normalisation
if (normalisation.method == "SCT") {
    seurat_obj <- SCTransform(seurat_obj, assay = "RNA", new.assay.name = "SCT", vars.to.regress = c("percent.mt"), variable.features.n = num.features, verbose = FALSE)
    #print(VariableFeatures(seurat_obj))
    #seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = num.features)
    #VariableFeatures(seurat_obj) <- features 
} else if (normalisation.method == "LogNormalize") {
    seurat_obj <- NormalizeData(object = seurat_obj, normalization.method = "LogNormalize", scale.factor = totalexpr)
    # Detection of variable genes across the single cells
    # FindVariableFeatures function identifies features that are outliers on a 'mean variability plot'.
    seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = num.features)

    ## Scaling:
    seurat_obj <- ScaleData(object = seurat_obj, vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = FALSE)
    print(seurat_obj)
}

# Open the pdf file for plotting
pdf(file = "Dispersion_plot.pdf", width = 13, height = 7)

#print(DefaultAssay(seurat_obj))

# Dispersion plot:
# Identify the 10 most highly variable genes
#top10 <- head(VariableFeatures(seurat_obj), 10)
# Plot variable features with and without labels
#dim(VariableFeatures(seurat_obj))
#plot1 <- VariableFeaturePlot(seurat_obj)
#plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
#CombinePlots(plots = list(plot1, plot2))

textplot(paste("\v \v Number of \n \v \v variable \n \v \v genes: \n \v \v", length(VariableFeatures(object = seurat_obj)), " \n  \n \v \v Number of \n \v \v cells: \n \v \v", length(colnames(x = seurat_obj))), halign = "center", valign = "center", cex = 2)

# Cell cycle stage scoring & PCA plot:
# Note: in the very beginning we read in the table and set s.genes and gm2.genes
# http://satijalab.org/seurat/cell_cycle_vignette.html#regress-out-cell-cycle-scores-during-data-scaling

# Check that there were some S or G2M genes in the list of variable genes:
if (length(s.genes[!is.na(match(s.genes, VariableFeatures(seurat_obj)))]) < 1 && length(g2m.genes[!is.na(match(g2m.genes, VariableFeatures(seurat_obj)))]) < 1) {
    # stop(paste('CHIPSTER-NOTE: ', "There were not enough cell cycle genes for correction in the list of variable genes."))
    # Write a log file (instead of ending with a Chipster-note, because we want the tool to finish and give the plot, when possible.)
    fileConn <- file("log.txt")
    writeLines(c("There are not enough cell cycle genes for correction in the list of variable genes."), fileConn)
    close(fileConn)
} else {
    # join layers before cell cycle scoring if global scaling normalization was used
    if (normalisation.method == "LogNormalize") {
        seurat_obj[["joined"]] <- JoinLayers(seurat_obj[["RNA"]])
        DefaultAssay(seurat_obj) <- "joined"
    }
    print(Layers(seurat_obj))
    print(DefaultAssay(seurat_obj))
    seurat_obj <- CellCycleScoring(seurat_obj,s.features = s.genes,g2m.features = g2m.genes, set.ident = TRUE)

    # Visualize in PCA:
    # PCA plot 1: without/before filtering cell cycle effect
    seurat_obj <- RunPCA(seurat_obj, features = c(s.genes, g2m.genes))
    plot1 <- DimPlot(seurat_obj) + ggtitle("PCA on cell cycle genes (no cell cycle regression)") # reduction = pca
    print("jejeje")
    # Cell cycle stage filtering:
    if (filter.cell.cycle != "no") {
        # Remove the cell cycle scores:

        # Option 1: remove all the difference:
        if (filter.cell.cycle == "all.diff") {
            # Scale again, but this time including also the cell cycle scores
            seurat_obj <- ScaleData(object = seurat_obj, vars.to.regress = c("S.Score", "G2M.Score", "nCount_RNA", "percent.mt"), verbose = FALSE)

            # PCA plot 2A: after filtering, all:
            seurat_obj <- RunPCA(object = seurat_obj, features = c(s.genes, g2m.genes), do.print = FALSE)
            ggtitle("After cell cycle correction (method: remove all)")
            plot2 <- DimPlot(seurat_obj) + ggtitle("After cell cycle correction (method: remove all)")
            CombinePlots(plots = list(plot1, plot2))
            # Option 2: regressing out the difference between the G2M and S phase scores:
        } else if (filter.cell.cycle == "diff.phases") {
            seurat_obj$CC.Difference <- seurat_obj$S.Score - seurat_obj$G2M.Score
            seurat_obj <- ScaleData(object = seurat_obj, vars.to.regress = c("CC.Difference", "nCount_RNA", "percent.mt"), verbose = FALSE)
            # PCA plot 2B: after filtering, difference:
            seurat_obj <- RunPCA(object = seurat_obj, features = c(s.genes, g2m.genes), do.print = FALSE)
            plot2 <- DimPlot(seurat_obj) + ggtitle("After cell cycle correction (method: G2M / S difference)")
            CombinePlots(plots = list(plot1, plot2))
        }
        # just plot the 1 PCA plot, if no filtering:
    } else {
        DimPlot(seurat_obj) # , plot.title = "PCA on cell cycle genes")
    }
    print(DefaultAssay(seurat_obj))
    # change back to original assay
    if (normalisation.method == "LogNormalize") {
        DefaultAssay(seurat_obj) <- "RNA"
    } else if (normalisation.method == "SCT") {
        DefaultAssay(seurat_obj) <- "SCT"
    }
}
dev.off() # close the pdf


# Save the Robj for the next tool
save(seurat_obj, file = "seurat_obj_combined.Robj")

## EOF
