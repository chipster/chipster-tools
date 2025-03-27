# started from comp-r-4-4-3 image

install.packages("BiocManager")
# required by RVAideMemoire (20 min?)
BiocManager::install("mixOmics")
# start with this, because it seems to have more difficult dependencies than other libraries
# requires apt package cmake
install.packages("RVAideMemoire")

install.packages("data.table")
install.packages("dplyr")
install.packages("ggplot2")
install.packages("vegan")
install.packages("plyr")
install.packages("ggpubr")

# 10 min?
BiocManager::install("microbiome")
BiocManager::install("DESeq2")

# something installed this already
# BiocManager::install("phyloseq")

install.packages("ashr")
