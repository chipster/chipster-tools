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

# CRAN packages and their dependencies
cranPackages = c(
		"amap", 
		"ape",
		"flexclust",
		"kohonen",
		"e1071",
		"fastICA",
		"aplpack",
		"corrgram",
		"deal",
		"outliers",
		"pvclust",
		"zoo",
		"scrime",
		"XML",
		"R2HTML",
		"moments",
		"snowfall",
		"sm",
		"rda",
		"flexmix", ## required by WECCA
		"MKmisc",
		"parallel",
		"doMC", # required by zinba
		"foreach", # required by zinba
		"quantreg", # required by zinba
		"R.utils", # required by zinba
		"tcltk",
		"ipred",
		"prodlim",
		"png",
		"gplots"
)

for (package in cranPackages) {
	smart.install.packages(package=package, mirror=repo.cran)
	library(package)
	detach("package:" + package, unload=TRUE)
}


# Bioconductor packages and their dependencies
bioconductorPackages = c(
		"RankProd",
		"Biobase",
		"IRanges",
		"AnnotationDbi",
		"affy",
		"affydata",
		"affyPLM",
		"affyQCReport",
		"annaffy",
		"annotate",
		"biomaRt",
		"DynDoc",
		"gcrma",
		"genefilter",
		"geneplotter",
		"hgu95av2.db",
		"limma",
		"marray",
		"multtest",
		"vsn",
		"xtable",
		"ctc",
		"oligo",
		"ssize",
		"LPE",
		"graph",
		"GlobalAncova",
		"impute",
		"idiogram",
		"GOstats",
		"beadarray",
		"simpleaffy",
		"globaltest",
		"locfit", # required by hdrcde 
		"ash", # required by hdrcde
		"ks", # required by hdrcde
		"hdrcde",
		"lumi",
		"prada",
		"siggenes",
		"plier",
		"beadarraySNP",
		"GEOquery",
		"ArrayExpress",
		"GeneCycle",
		"GeneNet",
		"GenomeGraphs",
		"MLInterfaces",
		"PAnnBuilder",
		"oligoClasses",
		"statmod",
		"vegan",
		"safe",
		"CGHbase",
		"CGHcall",
		"CGHregions",
		"DNAcopy",
		"MEDIPS",
		"microRNA",
		"RmiR",
		"ChIPpeakAnno",
		"rGADEM",
		"MotIV",
		"PICS",
		"seqLogo",
		"GenomicRanges",
		"Biostrings",
		"chipseq",
		"Rsamtools",
		"ShortRead",
		"lumiHumanIDMapping",
		"lumiMouseIDMapping",
		"lumiRatIDMapping",
		"edgeR",
		"rtracklayer",
		"maSigPro",
		"qvalue",
		"DESeq",
		"DESeq2",
		"DEXSeq",
		"RPA",
		"methylumi",
		"IlluminaHumanMethylation27k.db", # annotation package, not needed if all annotation packages from the repository are installed
		"VariantAnnotation",
		"TxDb.Hsapiens.UCSC.hg19.knownGene",
		"BSgenome.Hsapiens.UCSC.hg19",
		"rich",
		"BiodiversityR",
		"pegas",
		"labdsv",
		"sva",
		"Mfuzz",
		"WGCNA",
		"Heatplus",
		"RmiR.Hs.miRNA",
		"QDNAseq",
		
		# SNP 5.0 / 6.0
		"crlmm",
		"pd.mapping50k.hind240",
		"pd.mapping50k.xba240",
		"pd.mapping250k.nsp",
		"pd.mapping250k.sty",
		"pd.genomewidesnp.5",
		"pd.genomewidesnp.6",
		
		"genomewidesnp6Crlmm",
		"genomewidesnp5Crlmm",
		
		"human1mduov3bCrlmm",          #Illumina
		"human1mv1cCrlmm",             #Illumina
		"human370quadv3cCrlmm",        #Illumina
		"human370v1cCrlmm",            #Illumina
		"human550v3bCrlmm",            #Illumina
		"human610quadv1bCrlmm",        #Illumina
		"human650v3aCrlmm",            #Illumina
		"human660quadv1aCrlmm",        #Illumina
		"humanomni1quadv1bCrlmm",      #Illumina
		"humanomniexpress12v1bCrlmm"  #Illumina
)
#for (package in bioconductorPackages) {
#	smart.install.packages(bioconductor.package=package, mirror=repo.bioc)
#}



# Non-repo packages
nonRepoPackages = c(
		"http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R_libraries/NoWaves_0.4.tar.gz",
		"http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R_libraries/FruitFlyAgilent.db.tar.gz",
		"http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R_libraries/hgug4851a.db.tar.gz",
		"http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R_libraries/CGHtest_1.1.tar.gz",
		"http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R_libraries/CGHtestpar_0.0.tar.gz",
		"http://addictedtor.free.fr/packages/fpc_1.1-5.tar.gz",
		"http://www.math.utu.fi/projects/software/bio/ROTS_1.1.1.tar.gz",
		"http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R_libraries/A2R_0.0-4.tar.gz",
		"http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R_libraries/exon.pmcdf_1.1.tar.gz",
		"http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R_libraries/mouseexonpmcdf_1.1.tar.gz",
		"http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R_libraries/ratexonpmcdf_1.1.tar.gz",
		"http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R_libraries/intCNGEan_0.55.tar.gz",
		"http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R_libraries/CORNA_1.2.tar.gz",
		"http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R_libraries/WECCA_0.40.tar.gz",
		"http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R_libraries/zinba_2.02.03.tar.gz",
		"http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R_libraries/mgug4852a.db_1.0.0.tar.gz"
)

#for (package in nonRepoPackages) {
#	smart.install.packages(url.package=package)
#}




# Illumina annotation packages 
illuminaAnnotationPackages = c(
		"ftp://ctan.uib.no/pub/bioconductor/2.7/data/annotation/src/contrib/illuminaMousev1BeadID.db_1.8.0.tar.gz",
		"ftp://ctan.uib.no/pub/bioconductor/2.7/data/annotation/src/contrib/illuminaMousev2BeadID.db_1.8.0.tar.gz",
		"ftp://ctan.uib.no/pub/bioconductor/2.7/data/annotation/src/contrib/illuminaMousev1p1BeadID.db_1.8.0.tar.gz",
		"ftp://ctan.uib.no/pub/bioconductor/2.7/data/annotation/src/contrib/illuminaRatv1BeadID.db_1.8.0.tar.gz",
		"ftp://ctan.uib.no/pub/bioconductor/2.7/data/annotation/src/contrib/illuminaHumanv1BeadID.db_1.8.0.tar.gz",
		"ftp://ctan.uib.no/pub/bioconductor/2.7/data/annotation/src/contrib/illuminaHumanv2BeadID.db_1.8.0.tar.gz",
		"ftp://ctan.uib.no/pub/bioconductor/2.7/data/annotation/src/contrib/illuminaHumanv3BeadID.db_1.8.0.tar.gz",
		"ftp://ctan.uib.no/pub/bioconductor/2.7/data/annotation/src/contrib/illuminaHumanv4BeadID.db_1.8.0.tar.gz"
)

#for (package in illuminaAnnotationPackages) {
#	smart.install.packages(url.package=package);
#}


# Illumina annotation packages from bioconductor 
illuminaAnnotationBioconductorPackages = c(
		"illuminaHumanv1.db",
		"illuminaHumanv2.db",
		"illuminaHumanv3.db",
		"illuminaHumanv4.db",
		"illuminaMousev1.db",
		"illuminaMousev1p1.db",
		"illuminaMousev2.db",
		"illuminaRatv1.db",
		"org.Hs.eg.db",
		"org.Mm.eg.db",
		"org.Rn.eg.db",
		"org.Cf.eg.db",
		"TxDb.Hsapiens.UCSC.hg38.knownGene",
		"BSgenome.Hsapiens.UCSC.hg38",
		"BSgenome.Cfamiliaris.UCSC.canFam2",
		"BSgenome.Cfamiliaris.UCSC.canFam3",
		"PolyPhen.Hsapiens.dbSNP131"
)

#for (package in illuminaAnnotationBioconductorPackages) {
#	smart.install.packages(bioconductor.package=package, mirror=repo.bioc)
#}




# Install the whole annotation repository from Bioconductor
#smart.install.bioconductor.repo(repo.index = 3, mirror=repo.bioc) # for R 3.0.0, repo number 3 is annotations (might change)

# Install BrainArray custom CDF's. 
#smart.install.scavenge.web.packages("http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/latest/entrezg.asp", update=1)
#smart.install.scavenge.web.packages("http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/latest/tairg.asp", update=1)

# Install exon arrays from BrainArray. Parameter chiptype can be used to grep files with a certain substring. 
#smart.install.scavenge.web.packages("http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/latest/ense.asp", chiptype="ex", update=1)
#smart.install.scavenge.web.packages("http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/latest/ense.asp", chiptype="hta", update=1)

# Test if R.script has been linked to appropriate custom_cdf packages. Not perfect, but better than nothing.
# NOTE: below commands should match those listed above
#db.custom.packages <- smart.install.scavenge.web.packages("http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/latest/entrezg.asp", list.only=1)
#db.custom.packages <- c(db.custom.packages, smart.install.scavenge.web.packages("http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/latest/tairg.asp", list.only=1))
#db.custom.packages <- c(db.custom.packages, smart.install.scavenge.web.packages("http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/latest/ense.asp", chiptype="ex", list.only=1))
#db.custom.packages <- c(db.custom.packages, smart.install.scavenge.web.packages("http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/latest/ense.asp", chiptype="hta", list.only=1))
#check.affy.customnames(script.basename, db.custom.packages)
