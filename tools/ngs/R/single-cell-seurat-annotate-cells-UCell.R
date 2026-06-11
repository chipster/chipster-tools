# TOOL single-cell-seurat-annotate-cells-UCell.R: "Seurat v5 -Annotate cells with UCell" (You can use this tool to annotate the clusters.)
# INPUT seurat_obj.Robj: "Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL Plots.pdf
# PARAMETER OPTIONAL width: "Width of the output plots" TYPE INTEGER DEFAULT 10 (Width of the output plots in inches.)
# PARAMETER OPTIONAL height: "Height of the output plots" TYPE INTEGER DEFAULT 10 (Height of the output plots in inches.)
# PARAMETER  celltypes: "Cell types to plot" TYPE STRING DEFAULT "NK cells" (If you list multiple cell types, please use comma\(s\) \(,\) as a separator, e.g., \"L2/3 IT\,L4\".)
# PARAMETER  genesets: "Gene sets for celltypes" TYPE STRING (blahblah)
# PARAMETER OPTIONAL point.size: "Point size in tSNE and UMAP plots" TYPE DECIMAL DEFAULT 1 (Point size for the cluster plots.)
# PARAMETER OPTIONAL label.size: "Label size in the output plots" TYPE DECIMAL DEFAULT 4 (Label size in the output plots.)
# RUNTIME R-4.5.1-seurat5
# TOOLS_BIN ""


# 2026-04-06 JV


#Turn user given list into a vector

celltypes <- trimws(strsplit(celltypes, ",")[[1]])

genes <- trimws(strsplit(genesets, ",")[[1]])


library(UCell)
library(Seurat)
library(dplyr)
library(Matrix)
library(gplots)
library(ggplot2)


options(Seurat.object.assay.version = "v5")

## Load the R-Seurat-object (called seurat_obj)
load("seurat_obj.Robj")

if (exists("data.combined")) {
  seurat_obj <- data.combined
}

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



#Printtaa puuttuvat geenit, jos ei ole niin ota pois
print(missing)

seurat_obj <- AddModuleScore_UCell(seurat_obj, features=markers, name=NULL)


#Get the AddModuleScore -made columns:

old_names <- paste0(names(markers))
message(old_names)

score_mat <- seurat_obj@meta.data[, old_names]




#Get the max score for each cell (needed later because if its under 0 apparently, then its not the correct cell type???)
max_score <- apply(score_mat, 1, function(x) max(x))


#This line of code from https://stackoverflow.com/questions/17735859/for-each-row-return-the-column-name-of-the-largest-value
best_type <- colnames(score_mat)[apply(score_mat,1,which.max)]

#Remove the last number of cell type names created by AddModuleScore, e.g., "NK1" -> "NK"
best_type <- stringr::str_sub(best_type, start = 1, end = -1) 

seurat_obj$cell_type <- best_type


#Featureplot showing the score of each cell type (1 being max, 0 min)
#This is a sanity check for the researcher to see if the assigned cell type is about correct

pdf(file = "Plots.pdf", width = width, height = height)

p1 <- FeaturePlot(seurat_obj, label = T, label.size = label.size, reduction = "umap",  features = names(markers)[1:length(markers)])
p2 <- DimPlot(seurat_obj, group.by = "cell_type", label = T, label.size = label.size)

print(p1)
print(p2)

dev.off()



# EOF