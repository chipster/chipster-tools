
# CRAN packages and their dependencies
cranPackages = c(
		"Seurat",
		"dplyr",
		"pheatmap",
		"Matrix",
		"gplots",
		"cowplot",
		"ggplot2",
		"umap"
)

for (package in cranPackages) {
	install.packages(package=package)	
}

install.packages("BiocManager")
BiocManager::install("DESeq2")
BiocManager::install("MAST")
BiocManager::install("GEOquery")
BiocManager::install("scater")
BiocManager::install("mvoutlier")