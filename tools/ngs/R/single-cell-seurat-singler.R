# TOOL single-cell-seurat-singler.R: "SingleR cluster annotation" (Annotate your cell clusters using SingleR tool and CellDex annotation packages.)
# INPUT OPTIONAL seurat_obj.Robj: "Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL log.txt
# OUTPUT OPTIONAL singleR_annotation_plots.pdf
# OUTPUT OPTIONAL seurat_obj_singler_annotations.Robj
# OUTPUT OPTIONAL annotations_main.tsv
# OUTPUT OPTIONAL annotations_fine.tsv
# PARAMETER celldex.index: "CellDex reference to use" TYPE [HumanPrimaryCellAtlasData: "Human primary cell atlas", BlueprintEncodeData: "Blueprint ENCODE", MouseRNAseqData: "Mouse RNA-seq", ImmGenData: "Immunological Genome Project", DatabaseImmuneCellExpressionData: "Database of Immune Cell Expression", NovershternHematopoieticData: "Novershtern hematopoietic data", MonacoImmuneData: "Monaco Immune data"] DEFAULT MonacoImmuneData (Which CellDex reference to use for annotations.)
# PARAMETER OPTIONAL point.size: "Point size in plot" TYPE DECIMAL DEFAULT 1 (Point size for the cluster plots.)
# RUNTIME R-4.1.0-single-cell



# HumanPrimaryCellAtlasData, BlueprintEncodeData, MouseRNAseqData, ImmGenData, DatabaseImmuneCellExpressionData, NovershternHematopoieticData, MonacoImmuneData
# PARAMETER celldex.index "CellDex reference to use" TYPE [HumanPrimaryCellAtlasData: "Human primary cell atlas", BlueprintEncodeData: "Blueprint ENCODE", MouseRNAseqData: "Mouse RNA-seq", ImmGenData: "Immunological Genome Project", DatabaseImmuneCellExpressionData: "Database of Immune Cell Expression", NovershternHematopoieticData: "Novershtern hematopoietic data", MonacoImmuneData: "Monaco Immune data"] DEFAULT MonacoImmuneData (Which CellDex reference to use for annotations.)


# 2022-01-11 ML

# scRNAseq annotations with singleR
# Sources:
# 1. https://www.bioconductor.org/packages/release/bioc/vignettes/SingleR/inst/doc/SingleR.html
# 2. https://www.singlecellcourse.org/single-cell-rna-seq-analysis-using-seurat.html#cell-type-annotation-using-singler

library(Seurat)
library(dplyr)
library(Matrix)
library(gplots)
library(ggplot2)
library(celldex)
library(SingleR)

# Load the R-Seurat-object (called seurat_obj)
load("seurat_obj.Robj")

if (exists("data.combined")) {
  seurat_obj <- data.combined
}

# Let’s convert our Seurat object to single cell experiment (SCE) for convenience
sce <- as.SingleCellExperiment(DietSeurat(seurat_obj))


# HumanPrimaryCellAtlasData, BlueprintEncodeData, MouseRNAseqData, ImmGenData, DatabaseImmuneCellExpressionData, NovershternHematopoieticData, MonacoImmuneData
# HumanPrimaryCellAtlasData: "Human primary cell atlas", BlueprintEncodeData: "Blueprint/ENCODE", MouseRNAseqData: "Mouse RNA-seq", ImmGenData: "Immunological Genome Project", DatabaseImmuneCellExpressionData: "Database of Immune Cell Expression", NovershternHematopoieticData: "Novershtern hematopoietic data", MonacoImmuneData: "Monaco Immune data"


# annotation.data parameter = which annotation data, e.g. MonacoImmuneData
command <- paste("celldex::", celldex.index, "()", sep="")
ref <- eval(parse(text=command))


# These steps can take some time, < 1min
annotation.main <- SingleR(test = sce, assay.type.test = 1, ref = ref, labels = ref$label.main)
annotation.fine <- SingleR(test = sce, assay.type.test = 1, ref = ref, labels = ref$label.fine)

#write.table(annotations.main, file = "annotations_main.tsv", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
#write.table(annotations.fine, file = "annotations_fine.tsv", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)

write.table(table(annotation.main$pruned.labels), file = "annotations_main.tsv", sep = "\t", row.names = T, col.names = T, quote = F)
write.table(table(annotation.fine$pruned.labels), file = "annotations_fine.tsv", sep = "\t", row.names = T, col.names = T, quote = F)



# Back to Seurat obj
seurat_obj@meta.data$annotation.main <- annotation.main$pruned.labels
seurat_obj@meta.data$annotation.fine <- annotation.fine$pruned.labels


# Visualise, in pdf plots:
pdf(file = "singleR_annotation_plots.pdf",, width = 13, height = 7)
# Visualise:
seurat_obj <- SetIdent(seurat_obj, value = "annotation.main")
DimPlot(seurat_obj, label = T , repel = T, label.size = 3, pt.size = point.size) + NoLegend() + ggtitle('Main annotations')

seurat_obj <- SetIdent(seurat_obj, value = "annotation.fine")
DimPlot(seurat_obj, label = T , repel = T, label.size = 3, pt.size = point.size) + NoLegend() + ggtitle('Fine annotations')

plotScoreHeatmap(annotation.main)
plotDeltaDistribution(annotation.main, ncol = 3)
 
plotScoreHeatmap(annotation.fine)
plotDeltaDistribution(annotation.fine, ncol = 3)
dev.off() # close the pdf

# Save the Robj for the next tool
save(seurat_obj, file = "seurat_obj_singler_annotations.Robj")

## EOF
