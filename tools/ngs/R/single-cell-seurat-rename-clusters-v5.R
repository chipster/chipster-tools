# TOOL single-cell-seurat-rename-clusters-v5.R: "Seurat v5 -Rename clusters" (You can use this tool to rename the clusters.)
# INPUT seurat_obj.Robj: "Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL DimPlot.pdf
# PARAMETER OPTIONAL cluster_names.tsv: "Cluster name table in tsv format" TYPE STRING ()
# PARAMETER OPTIONAL celltypes: "Cell types to plot" TYPE STRING DEFAULT "NK cells" (If you list multiple cell types, please use comma\(s\) \(,\) as a separator, e.g., \"L2/3 IT\,L4\".)
# PARAMETER OPTIONAL genesets: "Gene sets for celltypes" TYPE STRING (blahblah)
# PARAMETER OPTIONAL point.size: "Point size in tSNE and UMAP plots" TYPE DECIMAL DEFAULT 1 (Point size for the cluster plots.)
# RUNTIME R-4.3.2-single-cell
# TOOLS_BIN ""


# 2021-12-27 ML
# 2023-12-15 IH
# 2026-04-06 JV

celltypes <- trimws(strsplit(celltypes, ",")[[1]])

genes <- trimws(strsplit(genesets, ",")[[1]])

library(Seurat)
library(dplyr)
library(Matrix)
library(gplots)
library(ggplot2)
options(Seurat.object.assay.version = "v5")

## OUTPUT OPTIONAL log.txt
## OUTPUT OPTIONAL seurat_obj_renamed.Robj
## OUTPUT OPTIONAL clusterPlotRenamed.pdf
## Load the R-Seurat-object (called seurat_obj)
load("seurat_obj.Robj")

if (exists("data.combined")) {
  seurat_obj <- data.combined
}

# Load the tsv with new cluster names:

if (exists("cluster_names")) {
clusternames <- read.table("cluster_names.tsv", sep = "\t", header = T, row.names = 1)

# Get the names from the second column:
new.cluster.ids <- clusternames[, 1]

# Some checks:
if (is.null(seurat_obj@commands$FindClusters)) {
  stop("CHIPSTER-NOTE: No cluster information in the Seurat object! Make sure you select an object that has gone through either Clustering or Integrated analysis of multiple samples tool.")
}
if (length(new.cluster.ids) != length(levels(seurat_obj))) {
  stop("CHIPSTER-NOTE: You need to give as input as many cluster names as there are clusters.")
}
if (!identical(row.names(clusternames), levels(seurat_obj))) {
  stop("CHIPSTER-NOTE: The cluster names = numbers in the input table need to be the same as in the Seurat object.")
}

# new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",  "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(seurat_obj)
seurat_obj <- RenameIdents(seurat_obj, new.cluster.ids)

} else {
seurat_obj$cell_type <- ""


#Get a gene and a group for it

genes_list <- c()
factors <- c()
nums <- 0

for(i in 1:length(genes)){
  if (genes[i] == "cut") {
    nums <- nums + 1
  } else {
    genes_list <- c(genes_list, genes[i])
    factors <- c(factors, nums)
  }
}


#x = pelkät geenin nimet, ilman seppiä
#f = toimii ryhmittäjänä nyt.


gene_groups <- split(x = genes_list, f = factors)
gene_groups

print("Rownames")
#message(rownames(seurat_obj))

print("Genes list")
#Features(seurat_obj[["RNA"]])

#Yhdistää geeniryhmät oikeisiin solutyyppeihin (Olettaa, että käyttäjä laitta oikeassa järjestyksessä)

markers <- setNames(gene_groups, celltypes)


markers
#Geeni pitää löytyä rivinimistä, muuten ei pysty laskemaan modulescorea -> Error


#Setdiff etsii ne alkiot, jotka ovat x:ssa mutteivat y:ssa
#Eli mitkä on markereissa ja mitkä eivät löydy sitten itse scRNA datasta
missing <- setdiff(
  unlist(markers, use.names = F),
  rownames(seurat_obj)
)

if (length(missing) > 0) {
  message("These genes are missing, please remove them:")
  message(paste0(missing, collapse = ", "))
  stop("CHIPSTER-NOTE: Some genes are missing")
}

#Printtaa puuttuvat geenit, jos ei ole niin ota pois
print(missing)

seurat_obj <- AddModuleScore(seurat_obj, features = markers, name = names(markers))

#Get the AddModuleScore -made columns:

old_names <- paste0(names(markers), seq_along(markers))
old_names

score_mat <- seurat_obj@meta.data[, old_names]


#Get the max score for each cell (needed later because if its under 0 apparently, then its not the correct cell type???)
max_score <- apply(score_mat, 1, function(x) max(x))


#This line of code from https://stackoverflow.com/questions/17735859/for-each-row-return-the-column-name-of-the-largest-value
best_type <- colnames(score_mat)[apply(score_mat,1,which.max)]

#Remove the last number of cell type names created by AddModuleScore, e.g., "NK1" -> "NK"
best_type <- stringr::str_sub(best_type, start = 1, end = -2) 



#Check if the addmodulescore is over 0 if so, then label with the highest score ("The most correct cell type"), and if its under 0, then unknown cell type
seurat_obj$cell_type <- ifelse(max_score > 0, best_type, "Other/Unknwon cell type")

pdf(file = "DimPlot.pdf", width = 13, height = 7)

print(DimPlot(seurat_obj, group.by = "cell_type"))

dev.off()

}


##Plot with renamed clusters:
#pdf(file = "clusterPlotRenamed.pdf")
#DimPlot(seurat_obj, reduction = "umap", label = TRUE, pt.size = point.size) + NoLegend()
#DimPlot(seurat_obj, reduction = "tsne", label = TRUE, pt.size = point.size)
# Number of cells in each cluster:
#cell_counts <- table(Idents(seurat_obj), seurat_obj$orig.ident)
#textplot(cell_counts, halign = "center", valign = "center", cex = 1)
#title(paste("Total number of cells: ", length(colnames(x = seurat_obj)), "\n Number of cells in each cluster:"))

#dev.off() # close the pdf

## Save the Robj for the next tool
##save(seurat_obj, file = "seurat_obj_renamed.Robj")

##



# EOF
