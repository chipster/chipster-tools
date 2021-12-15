
# CRAN packages and their dependencies

# this is only a copy of R-3.6.1-single-cell. The notes from the actual installation
# seem to be lost

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