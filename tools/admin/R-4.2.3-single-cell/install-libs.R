# CRAN packages and their dependencies

# these are needed to install some packages
install.packages("BiocManager")
install.packages("remotes")

# slow to install
install.packages("Seurat")
# 25 min
install.packages("Rfast2")
# 3 min
install.packages("Matrix")
# 11 min
BiocManager::install("DESeq2")
# 1 min
remotes::install_github("renozao/xbioc")

# install from pull request, because fix for NA values has not bee released yet
# < 4 min
remotes::install_github("meichendong/SCDC", ref = remotes::github_pull("31"))

# < 1 min
remotes::install_github("r-lib/systemfonts")
# requires r-lib/systemfonts
BiocManager::install("scater")

# < 1 min
BiocManager::install("multtest")
# requires multtest
# 3 min
install.packages("metap")

# fast to install
install.packages("dplyr")

install.packages("pheatmap")

install.packages("gplots")

install.packages("cowplot")

install.packages("ggplot2")
# > 3 min
install.packages("umap")

install.packages("tidyr")

install.packages("hdf5r")

# 1 min
BiocManager::install("MAST")
# 4 min
BiocManager::install("GEOquery")

BiocManager::install("mvoutlier")
# 2 min
BiocManager::install("celldex")
# 1 min
BiocManager::install("SingleR")
# 7 min
BiocManager::install("glmGamPoi")

# 2 min. This was replaced by Cellbender and shouldn't be installed next time
install.packages('SoupX')

install.packages("devtools")
# apt install git
devtools::install_github("cellgeni/schard")

install.packages("ggrepel")