
# CRAN packages and their dependencies

# this is only a copy of R-3.6.1-single-cell. The notes from the actual installation
# seem to be lost

cranPackages = c(
	"factoextra",
	"ggrepel"
)

for (package in cranPackages) {
	install.packages(package=package)	
}
