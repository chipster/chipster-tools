# TOOL single-cell-seurat-gene-plots.R: "Seurat v4 -Visualize genes" (Visualize for example selected cluster marker genes with violin and feature plot.)
# INPUT seurat_obj.Robj: "Seurat object" TYPE GENERIC
# INPUT OPTIONAL genes.txt: "Optional text file of the gene name(s)" TYPE GENERIC (The gene names\(s\) you wish to plot can also be given in the form of a text file, separated by comma. In case the text file is provided, the gene parameter is ignored.)
# OUTPUT OPTIONAL log.txt
# OUTPUT OPTIONAL seurat_obj_2.Robj
# OUTPUT OPTIONAL biomarker_plot.pdf
# OUTPUT OPTIONAL average_expressions.tsv 
# OUTPUT OPTIONAL percentage_of_cells_expressing.tsv
# PARAMETER OPTIONAL biomarker: "Gene name\(s\)" TYPE STRING DEFAULT "MS4A1, LYZ" (Name\(s\) of the biomarker gene to plot. If you list multiple gene names, use comma \(,\) as separator.)
# PARAMETER OPTIONAL normalisation.method: "Normalisation method used previously" TYPE [LogNormalize:"Global scaling normalization", SCT:"SCTransform"] DEFAULT LogNormalize (Which normalisation method was used in preprocessing, Global scaling normalization \(default, NormalizeData function used\) or SCTransform.)
# PARAMETER OPTIONAL point.size: "Point size in cluster plot" TYPE DECIMAL DEFAULT 1 (Point size for tSNE and UMAP plots.)
# PARAMETER OPTIONAL add.labels: "Add labels on top of clusters in plot" TYPE [TRUE: yes, FALSE: no] DEFAULT FALSE (Add cluster number on top of the cluster in UMAP plot.)
# PARAMETER OPTIONAL reduction.method: "Visualisation with tSNE, UMAP or PCA" TYPE [umap:UMAP, tsne:tSNE, pca:PCA] DEFAULT umap (Which dimensionality reduction plot to use.)
# PARAMETER OPTIONAL plotting.order.used: "Plotting order of cells based on expression" TYPE [TRUE:yes, FALSE:no] DEFAULT FALSE (Plot cells in the the order of expression. Can be useful to turn this on if cells expressing given feature are getting buried.)
# PARAMETER OPTIONAL output_aver_expr: "For each gene, list the average expression and percentage of cells expressing it in each cluster" TYPE [T: yes, F: no] DEFAULT F (Returns two tables: average expression and percentage of cells expressing the user defined genes in each cluster.)
# RUNTIME R-4.2.3-single-cell
# TOOLS_BIN ""



# 13.06.2017 ML
# 2018-01-11 ML update Seurat version to 2.2.0
# 2019-06-28 EK Add point size parameter for tSNE plot
# 2019-06-13 ML Seurat v3
# 2020-06-18 ML Add ridge plot
# 2021-10-04 ML Update to Seurat v4

# for UMAP:
library(reticulate)
Sys.setenv(RETICULATE_PYTHON = "/opt/chipster/tools-bin/miniconda3/envs/chipster_tools/bin/python")
#use_python("/opt/chipster/tools/miniconda3/envs/chipster_tools/bin/python")


library(Seurat)
library(dplyr)
library(Matrix)
library(gplots)
library(readr)  

# for the fileOk function
source(file.path(chipster.common.path,"tool-utils.R"))

# Load the R-Seurat-object (called seurat_obj)
load("seurat_obj.Robj")

if (exists("data.combined") ){
	seurat_obj <- data.combined
}


# In case some other type of assay is set:
if (normalisation.method == "SCT"){
  DefaultAssay(seurat_obj) <- "SCT"
}else{
 DefaultAssay(seurat_obj) <- "RNA"
}

# Use genes text file if provided, else the gene parameter is used
if (fileOk("genes.txt",0)) {
  genes = read_file("genes.txt")
  biomarker <- trimws(unlist(strsplit(genes,",")))
} else {
  biomarker <- trimws(unlist(strsplit(biomarker, ",")))
}

# Sanity check: are all of the requested genes available in the data (one missing allowed)
all.genes <- rownames(x = seurat_obj)
match(biomarker, all.genes)
# if more than one of the genes is not in the list, print error message:
if (sum(is.na((match(biomarker, all.genes)))) > 1) {  
  not.found <- (biomarker[is.na(match(biomarker, all.genes))==TRUE])
  not.found <- paste(not.found,collapse=",")
  stop(paste('CHIPSTER-NOTE: ', "The genes you requested were not found in this dataset:", not.found))
  }

# continue even if one gene is missing
if (!all(!is.na(match(biomarker, all.genes)))) { 
  not.found <- biomarker[is.na(match(biomarker, all.genes))==TRUE]
  print(paste("Continuing the visualization without the one gene not found: ", not.found))
  biomarker <- biomarker[!is.na(match(biomarker, all.genes))]
}

# open pdf
pdf(file="biomarker_plot.pdf", width=12, height=12) 

# Violin plot:
VlnPlot(seurat_obj, features = biomarker)

# Feature plot:
FeaturePlot(seurat_obj, features = biomarker, pt.size=point.size, reduction=reduction.method, label=add.labels, order=as.logical(plotting.order.used)) 	

# Ridge plot:
RidgePlot(seurat_obj, features = biomarker, ncol = 2)

# close the pdf
dev.off() 

# Average expression table
# If requested, return expression for an 'average' single cell in each cluster.
if (output_aver_expr == "T") {

  library(tidyr)

  a <- DotPlot(object = seurat_obj, features = biomarker)
  b <- a$data

  # Percentages:
  percentages <- select(b, features.plot, pct.exp, id)  %>% spread(id, pct.exp)
  row.names(percentages) <- percentages[,1] # pick the gene names as row names
  percentages <- round(percentages[,-1], digits=2) # remove gene name columns and round the numbers
  # Write to table
  write.table(percentages, file = "percentage_of_cells_expressing.tsv", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)

  # Average expressions:
  ave.expressions <- select(b, features.plot, avg.exp, id) %>% spread(id, avg.exp)
  row.names(ave.expressions) <- ave.expressions[,1] # pick the gene names as row names
  ave.expressions <- round(ave.expressions[,-1], digits=3) # remove gene name columns and round the numbers
  # Write to table
  write.table(ave.expressions, file = "average_expressions.tsv", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)


}

# EOF
