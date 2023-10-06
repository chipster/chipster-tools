# wget https://cran.r-project.org/src/base/R-4/R-4.3.1.tar.gz
# tar -xzf R-4.3.1.tar.gz
# cd R-4.3.1
# ./configure --with-x=no --with-pcre1 --prefix=/opt/chipster/tools/R-4.3.1-single-cell
# make
# make install

# CRAN packages and their dependencies

# these are needed to install some packages
install.packages("BiocManager")
install.packages("remotes")

# slow to install
# > 54 min
# install from github, because seurat5 is still in beta
remotes::install_github("satijalab/seurat", "seurat5", quiet = TRUE)
#install.packages("Seurat")
# 15:58 - 
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
