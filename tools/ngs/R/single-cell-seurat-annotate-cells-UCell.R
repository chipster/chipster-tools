# TOOL single-cell-seurat-annotate-cells-UCell.R: "Seurat v5 - BETA Annotate cells with UCell" (You can use this tool to annotate the clusters with own cell types and gene sets.)
# INPUT seurat_obj.Robj: "Seurat object. Has to be pre-processed so that it contains UMAP information" TYPE GENERIC
# INPUT OPTIONAL celltypes_markers.tsv: "Cell types and gene sets" TYPE GENERIC
# OUTPUT UCell_Plots.pdf
# OUTPUT seurat_obj_annotated_UCell.Robj
# PARAMETER OPTIONAL celltypes: "Cell types to plot" TYPE STRING DEFAULT "T cells, B cells, NK cells, Monocytes" (Please use comma\(s\) \(,\) as a separator, e.g., \T Cells\, B cells\. Minimum of 2 cell types is required because UCell is based on Mann-Whitney U test) 
# PARAMETER OPTIONAL genesets: "Gene sets for cell types" TYPE STRING DEFAULT "CD3D, CD3E, IL7R, sep, MS4A1, CD79A, CD79B, sep, NKG7, GNLY, PRF1, sep, CD14, LST1, S100A8" (Gene sets for cell types in the same order as celltypes. For example if you have T cells, B cells, NK cells, Monocytes, please input first the geneset for T cells, then use word "sep" and then type in geneset for B cells and so on. If you list multiple gene sets, please use comma\(s\) \(,\) as a separator, e.g., \CD3D\, CD3E\, IL7R\, sep\, MS4A1\, CD79A\, CD79B\, sep\, NKG7\, GNLY\, PRF1\, sep\, CD14\, LST1\, S100A8\. ) 
# PARAMETER OPTIONAL width: "Width of the output plots" TYPE INTEGER DEFAULT 10 (Width of the output plots in inches.)
# PARAMETER OPTIONAL height: "Height of the output plots" TYPE INTEGER DEFAULT 10 (Height of the output plots in inches.)
# PARAMETER OPTIONAL point.size: "Point size in tSNE and UMAP plots" TYPE DECIMAL DEFAULT 1 (Point size for the cluster plots.)
# PARAMETER OPTIONAL label.size: "Label size in the output plots" TYPE DECIMAL DEFAULT 4 (Label size for cluster numbers or cell type names on top of UMAP. If you don't want any labels, set this to 0.)
# RUNTIME R-4.5.1-seurat5
# TOOLS_BIN ""

chipster.threads.max <- as.numeric(chipster.threads.max)

# 2026-04-06 JV 

#Load needed packages
library(UCell)
library(Seurat)
library(dplyr)
library(Matrix)
library(gplots)
library(ggplot2)


options(Seurat.object.assay.version = "v5")


# Try catch checking if the correct file has been inputted

try_catch_load <- function(file) {
  tryCatch({
    suppressWarnings(load(file, envir = globalenv()))
  }, error = function(e) {
    stop("CHIPSTER-NOTE: Failed to load the file, try swapping the files", call. = FALSE)
  })
}

# Load the R-Seurat-object (called seurat_obj)

try_catch_load("seurat_obj.Robj")

# Make a new column called cell_type and fill it with "" 
seurat_obj$cell_type <- ""

if (file.exists("celltypes_markers.tsv")) {
  file <- read.table("celltypes_markers.tsv", sep = "\t", header = T)

  celltypes <- file[,1]
  genes <- file[,2]

  markers <- split(x = genes, f = celltypes)

  # Make the gene list in the proper format
  markers <- lapply(markers, function(x) {
    trimws(strsplit(x, ",")[[1]])
  })

} else {

celltypes <- trimws(strsplit(celltypes, ",")[[1]])
genes <- trimws(strsplit(genesets, ",")[[1]])


# Empty vectors for the for loop
genes_list <- c()
factors <- c()
nums <- 0

# This basically assigns the cell type and its marker-genes so that the genelist gets separated by "sep"
for(i in 1:length(genes)){
  if (genes[i] == "sep") {
    nums <- nums + 1
  } else {
    genes_list <- c(genes_list, genes[i])
    factors <- c(factors, nums)
  }
}

# This groups the genes based on their celltype group

gene_groups <- split(x = genes_list, f = factors)

# This combines the gene groups with the correct cell types (Assumes that the user has put them in the correct order)

markers <- setNames(gene_groups, celltypes)

}

# The gene has to be found in the rownames, otherwise it can't calculate the modulescore -> Error

#Find genes not found in the seurat object rownames

gene_names <- unlist(markers, use.names = F)
gene_names <- sub("[+-]$", "", gene_names)

missing <- setdiff(
  gene_names,
  rownames(seurat_obj)
)


if (length(missing) > 0) {
  message("The following genes were not found in the Seurat object. Expression will be imputed as 0. Genes not found: ", paste(missing, collapse = ", "))
}


# Calculate UCell module score for each cell type and add it to the Seurat object. 
seurat_obj <- AddModuleScore_UCell(seurat_obj, features=markers, name=NULL, ncore = chipster.threads.max)


# Get scores out of the seurat object
score_mat <- seurat_obj@meta.data[, names(markers)]


# This line of code from https://stackoverflow.com/questions/17735859/for-each-row-return-the-column-name-of-the-largest-value
# Gets the cell type that has the highest module score

best_type <- colnames(score_mat)[apply(score_mat,1,which.max)]

seurat_obj$cell_type <- best_type



# Assign a cell type to each cluster based on the majority of cells in that cluster


seurat_table <- table(seurat_obj$seurat_clusters, seurat_obj$cell_type)

type <- apply(seurat_table, 1, function(x) names(which.max(x)))


# Check that cluster numbers on the seurat object and on the table actually match


match <- all(seurat_obj$seurat_clusters == names(type[as.character(seurat_obj$seurat_clusters)]))

if (!match) {
  "CHIPSER-NOTE: Cluster number on the seurat object does not match with the seurat_table. Wrong cell types would be annotated. Stopping process..."
} else {
  print("Cluster numbers match, assigning a cell type per cluster")
}

# head(type, 20)
seurat_obj$cluster_celltype <- as.vector(type[as.character(seurat_obj$seurat_clusters)])



# Featureplot showing the score of each cell type (1 being max, 0 min)
# This is a sanity check for the researcher to see if the assigned cell type is about correct

pdf(file = "UCell_Plots.pdf", width = width, height = height)

p1 <- FeaturePlot(seurat_obj, label = T, label.size = label.size, pt.size = point.size, reduction = "umap",  features = names(markers)[1:length(markers)])+
  labs(color = "Module scores")


p2 <- DimPlot(seurat_obj, group.by = "cell_type", label = T, pt.size = point.size, label.size = label.size)+
  ggtitle("Assigned cell types based on highest UCell module score")+
  labs(color = "Cell types")

p3 <- DimPlot(seurat_obj, group.by = "cluster_celltype", label = T, pt.size = point.size, label.size = label.size)+
  ggtitle("Assigned cell types based on majority of cells in cluster")+
  labs(color = "Cell types")


print(p1)
print(p2)
print(p3)

dev.off()

save(seurat_obj, file = "seurat_obj_annotated_UCell.Robj")

# EOF