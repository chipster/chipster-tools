
# CRAN packages and their dependencies
install.packages("Seurat")
install.packages("dplyr")
install.packages("pheatmap")
install.packages("Matrix")
install.packages("gplots")
install.packages("cowplot")
install.packages("ggplot2")
install.packages("umap")

install.packages("BiocManager")
BiocManager::install("DESeq2")
BiocManager::install("MAST")
BiocManager::install("GEOquery")
BiocManager::install("scater")
BiocManager::install("mvoutlier")
BiocManager::install('multtest')

# requires multtest
install.packages('metap')