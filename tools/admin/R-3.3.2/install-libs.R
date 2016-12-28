# This install script is run in bare R installation and 
# it installs all packages required to run Chipster.
# The script uses install functions that check each package
# before installation, meaning that it can be rerun if needed
# and only missing packages are installed.

# Determine the path of the executing script
initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)
# Make path name absolute
script.basename <- normalizePath(script.basename)

# Use smart.* install utility functions
# They skip all packages that already have been installed
source(paste(script.basename, "/smip.R", sep=""));

# Configure paths and repos (change if you need)
repo.cran <- "http://ftp.acc.umu.se/mirror/CRAN/"
repo.bioc <- "http://www.bioconductor.org"

#check where this script resides
#relative.script.dir <- dirname(parent.frame(2)$ofile)
#absolute.script.dir <- normalizePath(relative.script.dir)
#source(paste(absolute.script.dir, "/smip.R", sep=""))

# Use smart.* install utility functions
# They skip all packages that already have been installed
#source("smip.R")

#Unlike R, RScript does not seem to load the method-package, why some try-catches can crash
library(methods)

# Bioconductor packages and their dependencies
bioconductorPackages = c(
		"Biobase",
		"QDNAseq"
)
		
for (package in bioconductorPackages) {
	smart.install.packages(bioconductor.package=package, mirror=repo.bioc)
	#library(package, character.only = TRUE)
	#detach(paste("package:", package, sep = ""), character.only = TRUE, unload=TRUE)
}


