# wget https://cran.r-project.org/src/base/R-4/R-4.3.2.tar.gz
# tar -xzf R-4.3.2.tar.gz 
# cd R-4.3.2
# ./configure --with-x=no --with-pcre1 --prefix=/opt/chipster/tools/R-4.3.2-single-cell
# make
# make install

# First, install all apt package in chipster-openshift/kustomize/builds/comp-r-4-3-2-single-cell!

# you can run R commands from bash to track time usage, but you have to provide the repo for install.packages() 
# time ./R -e 'install.packages("Seurat", repos = "https://ftp.acc.umu.se/mirror/CRAN/")'

# these are needed to install some packages
install.packages("BiocManager")
install.packages("remotes")

# slow to install
# 55 min
install.packages('Seurat')

# "Seurat does not require, but makes use of, packages developed by other labs that can substantially enhance speed and performance"
setRepositories(ind = 1:3, addURLs = c('https://satijalab.r-universe.dev', 'https://bnprks.r-universe.dev/'))
# 40 min?
# BPCells requires "apt-get install -y libhdf5-dev"
install.packages(c("BPCells", "presto", "glmGamPoi"))

# packages which required manually installed dependencies

# < 1 min
remotes::install_github("r-lib/systemfonts")
# requires r-lib/systemfonts
BiocManager::install("scater")

# < 1 min
BiocManager::install("multtest")
# requires multtest
# 3 min
install.packages("metap")


# slow packages from install.packages

# 29 min
install.packages("Rfast2")
# 3 min
install.packages("Matrix")
# 1 min
install.packages("umap")
# 2 min. This was replaced by Cellbender and shouldn't be installed next time
install.packages('SoupX')


# BiocManager

# 1 min
BiocManager::install("DESeq2")
# something installed this already
#BiocManager::install("glmGamPoi")
# 6 min
BiocManager::install("GEOquery")
# 6 min
BiocManager::install("celldex")
# 1 min
BiocManager::install("MAST")
# 2 min
BiocManager::install("SingleR")
# 17 s
BiocManager::install("mvoutlier")


# remotes

# install from pull request, because fix for NA values has not bee released yet
# 1 min
remotes::install_github("meichendong/SCDC", ref = remotes::github_pull("31"))
# 1 min
remotes::install_github("renozao/xbioc")


# fast to install
install.packages("dplyr")

install.packages("pheatmap")

install.packages("gplots")

install.packages("cowplot")

install.packages("ggplot2")

install.packages("tidyr")

install.packages("hdf5r")





