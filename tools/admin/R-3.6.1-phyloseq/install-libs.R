# wget https://ftp.acc.umu.se/mirror/CRAN/src/base/R-3/R-3.6.1.tar.gz
# tar -xzf R-3.6.1.tar.gz
# cd R-3.6.1
# ./configure --prefix=/opt/chipster/tools/R-3.6.1-phyloseq
# # ./configure --prefix=/opt/chipster/tools/R-3.6.1-phyloseq --with-pcre1
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
# I just happened to install in this order
BiocManager::install('microbiome')
install.packages("plyr")

BiocManager::install("DESeq2")
BiocManager::install("phyloseq")