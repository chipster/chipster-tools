# The following commands were used to build R to tools-bin
# 
# cd /opt/chipster/tools/temp
# wget https://ftp.acc.umu.se/mirror/CRAN/src/base/R-4/R-4.1.1.tar.gz
# tar -xzf R-4.1.1.tar.gz
# cd R-4.1.1
# ./configure --with-x=no --with-pcre1 --prefix=/opt/chipster/tools/R-4.1.1
# make
# make install

# CRAN packages and their dependencies
cranPackages = c(
)

for (package in cranPackages) {
	install.packages(package=package)	
}
