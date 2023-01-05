# wget https://ftp.acc.umu.se/mirror/CRAN/src/base/R-4/R-4.2.0.tar.gz
# tar -xf R-4.2.0.tar.gz
# cd R-4.2.0
# ./configure  --with-x=no --with-pcre1 --prefix=/opt/chipster/tools/R-4.2.0-phyloseq
# make
# make install

install.packages("BiocManager")
# required by RVAideMemoire
BiocManager::install('mixOmics')
# start with this, because it seems to have more difficult dependencies than other libraries
install.packages("RVAideMemoire")

install.packages("data.table")
install.packages("dplyr")
install.packages("ggplot2")
install.packages("vegan")
install.packages("plyr")
install.packages("ggpubr")

BiocManager::install('microbiome')
BiocManager::install("DESeq2")

# something installed this already
# BiocManager::install("phyloseq")

install.packages("Rfast2")
