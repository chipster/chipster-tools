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
source(paste(script.basename, "/smip.R", sep = ""))
# Configure paths and repos (change if you need)
repo.cran <- "http://ftp.sunet.se/pub/lang/CRAN"
repo.bioc <- "http://www.bioconductor.org"

# check where this script resides
# relative.script.dir <- dirname(parent.frame(2)$ofile)
# absolute.script.dir <- normalizePath(relative.script.dir)
# source(paste(absolute.script.dir, "/smip.R", sep=""))

# Use smart.* install utility functions
# They skip all packages that already have been installed
# source("smip.R")

# Unlike R, RScript does not seem to load the method-package, why some try-catches can crash
library(methods)

# Install packages, and their dependencies, from CRAN
smart.install.packages(package = "amap", mirror = repo.cran)
smart.install.packages(package = "ape", mirror = repo.cran)
smart.install.packages(package = "flexclust", mirror = repo.cran)
smart.install.packages(package = "kohonen", mirror = repo.cran)
smart.install.packages(package = "e1071", mirror = repo.cran)
smart.install.packages(package = "fastICA", mirror = repo.cran)
smart.install.packages(package = "aplpack", mirror = repo.cran)
smart.install.packages(package = "corrgram", mirror = repo.cran)
smart.install.packages(package = "deal", mirror = repo.cran)
smart.install.packages(package = "outliers", mirror = repo.cran)
smart.install.packages(package = "pvclust", mirror = repo.cran)
smart.install.packages(package = "zoo", mirror = repo.cran)
smart.install.packages(package = "scrime", mirror = repo.cran)
smart.install.packages(package = "XML", mirror = repo.cran)
smart.install.packages(package = "R2HTML", mirror = repo.cran)
smart.install.packages(package = "moments", mirror = repo.cran)
smart.install.packages(package = "snowfall", mirror = repo.cran)
smart.install.packages(package = "sm", mirror = repo.cran)
smart.install.packages(package = "rda", mirror = repo.cran)
smart.install.packages(package = "flexmix", mirror = repo.cran) ## required by WECCA
smart.install.packages(package = "MKmisc", mirror = repo.cran)
smart.install.packages(package = "multicore", mirror = repo.cran) # required by zinba
smart.install.packages(package = "doMC", mirror = repo.cran) # required by zinba
smart.install.packages(package = "foreach", mirror = repo.cran) # required by zinba
smart.install.packages(package = "quantreg", mirror = repo.cran) # required by zinba
smart.install.packages(package = "R.utils", mirror = repo.cran) # required by zinba
smart.install.packages(package = "tcltk", mirror = repo.cran)
smart.install.packages(package = "ipred", mirror = repo.cran)
smart.install.packages(package = "prodlim", mirror = repo.cran)
smart.install.packages(package = "png", mirror = repo.cran)

# Install packages, and their dependencies, from Bioconductor
smart.install.packages(bioconductor.package = "RankProd", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "Biobase", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "IRanges", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "AnnotationDbi", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "affy", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "affydata", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "affyPLM", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "affyQCReport", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "annaffy", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "annotate", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "biomaRt", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "Biostrings", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "DynDoc", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "gcrma", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "genefilter", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "geneplotter", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "GenomicRanges", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "hgu95av2.db", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "limma", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "marray", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "multtest", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "vsn", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "xtable", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "ctc", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "oligo", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "ssize", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "LPE", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "graph", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "affyQCReport", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "GlobalAncova", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "impute", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "idiogram", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "GOstats", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "beadarray", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "simpleaffy", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "globaltest", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "geneplotter", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "biomaRt", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "locfit", mirror = repo.bioc) # required by hdrcde
smart.install.packages(bioconductor.package = "ash", mirror = repo.bioc) # required by hdrcde
smart.install.packages(bioconductor.package = "ks", mirror = repo.bioc) # required by hdrcde
smart.install.packages(bioconductor.package = "hdrcde", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "lumi", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "prada", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "siggenes", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "plier", mirror = repo.bioc)
# smart.install.packages(bioconductor.package="cosmo", mirror=repo.bioc) # missing
smart.install.packages(bioconductor.package = "beadarraySNP", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "GEOquery", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "ArrayExpress", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "GeneCycle", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "GeneNet", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "GenomeGraphs", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "MLInterfaces", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "GOstats", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "PAnnBuilder", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "oligoClasses", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "statmod", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "vegan", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "safe", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "CGHbase", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "CGHcall", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "CGHregions", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "DNAcopy", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "MLInterfaces", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "MEDIPS", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "microRNA", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "RmiR", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "biomaRt", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "ChIPpeakAnno", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "rGADEM", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "MotIV", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "PICS", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "seqLogo", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "GenomicRanges", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "Biostrings", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "ChIPpeakAnno", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "chipseq", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "Rsamtools", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "ShortRead", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "lumiHumanIDMapping", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "lumiMouseIDMapping", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "lumiRatIDMapping", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "edgeR", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "rtracklayer", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "maSigPro", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "qvalue", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "DESeq", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "DESeq2", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "DEXSeq", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "RPA", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "methylumi", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "IlluminaHumanMethylation450k.db", mirror = repo.bioc) # annotation package, not needed if all annotation packages from the repository are installed
smart.install.packages(bioconductor.package = "IlluminaHumanMethylation27k.db", mirror = repo.bioc) # annotation package, not needed if all annotation packages from the repository are installed
smart.install.packages(bioconductor.package = "VariantAnnotation", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "TxDb.Hsapiens.UCSC.hg19.knownGene", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "DEXSeq", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "BSgenome.Hsapiens.UCSC.hg19", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "rich", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "BiodiversityR", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "pegas", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "labdsv", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "sva", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "Mfuzz", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "WGCNA", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "Heatplus", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "RmiR.Hs.miRNA", mirror = repo.bioc)

# SNP 5.0 / 6.0
smart.install.packages(bioconductor.package = "crlmm", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "oligo", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "pd.mapping50k.hind240", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "pd.mapping50k.xba240", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "pd.mapping250k.nsp", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "pd.mapping250k.sty", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "pd.genomewidesnp.5", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "pd.genomewidesnp.6", mirror = repo.bioc)

smart.install.packages(bioconductor.package = "genomewidesnp6Crlmm", mirror = repo.bioc)
smart.install.packages(bioconductor.package = "genomewidesnp5Crlmm", mirror = repo.bioc)

smart.install.packages(bioconductor.package = "human1mduov3bCrlmm", mirror = repo.bioc) # Illumina
smart.install.packages(bioconductor.package = "human1mv1cCrlmm", mirror = repo.bioc) # Illumina
smart.install.packages(bioconductor.package = "human370quadv3cCrlmm", mirror = repo.bioc) # Illumina
smart.install.packages(bioconductor.package = "human370v1cCrlmm", mirror = repo.bioc) # Illumina
smart.install.packages(bioconductor.package = "human550v3bCrlmm", mirror = repo.bioc) # Illumina
smart.install.packages(bioconductor.package = "human610quadv1bCrlmm", mirror = repo.bioc) # Illumina
smart.install.packages(bioconductor.package = "human650v3aCrlmm", mirror = repo.bioc) # Illumina
smart.install.packages(bioconductor.package = "human660quadv1aCrlmm", mirror = repo.bioc) # Illumina
smart.install.packages(bioconductor.package = "humanomni1quadv1bCrlmm", mirror = repo.bioc) # Illumina
smart.install.packages(bioconductor.package = "humanomniexpress12v1bCrlmm", mirror = repo.bioc) # Illumina


# Install non-repo packages
smart.install.packages(url.package = "http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R_libraries/NoWaves_0.4.tar.gz")
smart.install.packages(url.package = "http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R_libraries/FruitFlyAgilent.db.tar.gz")
smart.install.packages(url.package = "http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R_libraries/hgug4851a.db.tar.gz")
smart.install.packages(url.package = "http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R_libraries/CGHtest_1.1.tar.gz")
smart.install.packages(url.package = "http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R_libraries/CGHtestpar_0.0.tar.gz")
smart.install.packages(url.package = "http://addictedtor.free.fr/packages/fpc_1.1-5.tar.gz")
smart.install.packages(url.package = "http://www.math.utu.fi/projects/software/bio/ROTS_1.1.1.tar.gz")
smart.install.packages(url.package = "http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R_libraries/A2R_0.0-4.tar.gz")
smart.install.packages(url.package = "http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R_libraries/exon.pmcdf_1.1.tar.gz")
smart.install.packages(url.package = "http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R_libraries/mouseexonpmcdf_1.1.tar.gz")
smart.install.packages(url.package = "http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R_libraries/ratexonpmcdf_1.1.tar.gz")
smart.install.packages(url.package = "http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R_libraries/intCNGEan_0.55.tar.gz")
smart.install.packages(url.package = "http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R_libraries/CORNA_1.2.tar.gz")
smart.install.packages(url.package = "http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R_libraries/WECCA_0.40.tar.gz")
smart.install.packages(url.package = "http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R_libraries/zinba_2.02.03.tar.gz")
smart.install.packages("QDNAseq", mirror = "http://www.bioconductor.org/packages/2.14/bioc", update = 1) # This package is in Bioconductor only since version 2.14 (which is for R-3.1), hence the "forced" install instead of regular Bioconductor installation (which for R-3.0 uses Bioconductor 2.13).
smart.install.packages(url.package = "http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R_libraries/mgug4852a.db_1.0.0.tar.gz")


# affy_20 for R-2.12"
# smart.install.packages(url.package="http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R_libraries/affy_20/aragene11stattairgcdf_17.0.0.tar.gz")
# smart.install.packages(url.package="http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R_libraries/affy_20/aragene11stattairgprobe_17.0.0.tar.gz")
# smart.install.packages(url.package="http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R_libraries/affy_20/aragene10stattairgcdf_17.0.0.tar.gz")
# smart.install.packages(url.package="http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R_libraries/affy_20/aragene10stattairgprobe_17.0.0.tar.gz")
# smart.install.packages(url.package="http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R_libraries/affy_20/hugene20sthsentrezgcdf_17.0.0.tar.gz")
# smart.install.packages(url.package="http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R_libraries/affy_20/hugene20sthsentrezgprobe_17.0.0.tar.gz")
# smart.install.packages(url.package="http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R_libraries/affy_20/hugene21sthsentrezgcdf_17.0.0.tar.gz")
# smart.install.packages(url.package="http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R_libraries/affy_20/hugene21sthsentrezgprobe_17.0.0.tar.gz")
# smart.install.packages(url.package="http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R_libraries/affy_20/aragene10statentrezgcdf_17.0.0.tar.gz")
# smart.install.packages(url.package="http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R_libraries/affy_20/aragene10statentrezgprobe_17.0.0.tar.gz")
# smart.install.packages(url.package="http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R_libraries/affy_20/aragene11statentrezgprobe_17.0.0.tar.gz")
# smart.install.packages(url.package="http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R_libraries/affy_20/aragene11statentrezgcdf_17.0.0.tar.gz")
# smart.install.packages(url.package="http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R_libraries/affy_20/ragene20strnentrezgcdf_17.0.0.tar.gz")
# smart.install.packages(url.package="http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R_libraries/affy_20/ragene20strnentrezgprobe_17.0.0.tar.gz")
# smart.install.packages(url.package="http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R_libraries/affy_20/ragene21strnentrezgprobe_17.0.0.tar.gz")
# smart.install.packages(url.package="http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R_libraries/affy_20/ragene21strnentrezgcdf_17.0.0.tar.gz")
# smart.install.packages(url.package="http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R_libraries/affy_20/mogene20stmmentrezgcdf_17.0.0.tar.gz")
# smart.install.packages(url.package="http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R_libraries/affy_20/mogene20stmmentrezgprobe_17.0.0.tar.gz")
# smart.install.packages(url.package="http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R_libraries/affy_20/mogene21stmmentrezgcdf_17.0.0.tar.gz")
# smart.install.packages(url.package="http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R_libraries/affy_20/mogene21stmmentrezgprobe_17.0.0.tar.gz")
# smart.install.packages(url.package="http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R_libraries/affy_20/zebgene10stdrentrezgcdf_17.0.0.tar.gz")
# smart.install.packages(url.package="http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R_libraries/affy_20/zebgene10stdrentrezgprobe_17.0.0.tar.gz")
# smart.install.packages(url.package="http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R_libraries/affy_20/zebgene11stdrentrezgprobe_17.0.0.tar.gz")
# smart.install.packages(url.package="http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/R_libraries/affy_20/zebgene11stdrentrezgcdf_17.0.0.tar.gz")

# Install Illumina annotation packages.
# smart.install.packages(url.package="ftp://ftp.au.freebsd.org/pub/bioconductor/packages/2.8/data/annotation/src/contrib/illuminaMousev1BeadID.db_1.8.0.tar.gz");
# smart.install.packages(url.package="ftp://ftp.au.freebsd.org/pub/bioconductor/packages/2.8/data/annotation/src/contrib/illuminaMousev2BeadID.db_1.8.0.tar.gz");
# smart.install.packages(url.package="ftp://ftp.au.freebsd.org/pub/bioconductor/packages/2.8/data/annotation/src/contrib/illuminaMousev1p1BeadID.db_1.8.0.tar.gz");
# smart.install.packages(url.package="ftp://ftp.au.freebsd.org/pub/bioconductor/packages/2.8/data/annotation/src/contrib/illuminaRatv1BeadID.db_1.8.0.tar.gz");
# smart.install.packages(url.package="ftp://ftp.au.freebsd.org/pub/bioconductor/packages/2.8/data/annotation/src/contrib/illuminaHumanv1BeadID.db_1.8.0.tar.gz");
# smart.install.packages(url.package="ftp://ftp.au.freebsd.org/pub/bioconductor/packages/2.8/data/annotation/src/contrib/illuminaHumanv2BeadID.db_1.8.0.tar.gz");
# smart.install.packages(url.package="ftp://ftp.au.freebsd.org/pub/bioconductor/packages/2.8/data/annotation/src/contrib/illuminaHumanv3BeadID.db_1.8.0.tar.gz");
# smart.install.packages(url.package="ftp://ftp.au.freebsd.org/pub/bioconductor/packages/2.8/data/annotation/src/contrib/illuminaHumanv4BeadID.db_1.8.0.tar.gz");
smart.install.packages(url.package = "ftp://ctan.uib.no/pub/bioconductor/2.7/data/annotation/src/contrib/illuminaMousev1BeadID.db_1.8.0.tar.gz")
smart.install.packages(url.package = "ftp://ctan.uib.no/pub/bioconductor/2.7/data/annotation/src/contrib/illuminaMousev2BeadID.db_1.8.0.tar.gz")
smart.install.packages(url.package = "ftp://ctan.uib.no/pub/bioconductor/2.7/data/annotation/src/contrib/illuminaMousev1p1BeadID.db_1.8.0.tar.gz")
smart.install.packages(url.package = "ftp://ctan.uib.no/pub/bioconductor/2.7/data/annotation/src/contrib/illuminaRatv1BeadID.db_1.8.0.tar.gz")
smart.install.packages(url.package = "ftp://ctan.uib.no/pub/bioconductor/2.7/data/annotation/src/contrib/illuminaHumanv1BeadID.db_1.8.0.tar.gz")
smart.install.packages(url.package = "ftp://ctan.uib.no/pub/bioconductor/2.7/data/annotation/src/contrib/illuminaHumanv2BeadID.db_1.8.0.tar.gz")
smart.install.packages(url.package = "ftp://ctan.uib.no/pub/bioconductor/2.7/data/annotation/src/contrib/illuminaHumanv3BeadID.db_1.8.0.tar.gz")
smart.install.packages(url.package = "ftp://ctan.uib.no/pub/bioconductor/2.7/data/annotation/src/contrib/illuminaHumanv4BeadID.db_1.8.0.tar.gz")
# Install the whole annotation repository from Bioconductor
smart.install.bioconductor.repo(repo.index = 3, mirror = repo.bioc) # for R 3.0.0, repo number 3 is annotations (might change)

# Install BrainArray custom CDF's.
smart.install.scavenge.web.packages("http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/latest/entrezg.asp", update = 1)
smart.install.scavenge.web.packages("http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/latest/tairg.asp", update = 1)
# Install exon arrays from BrainArray. Parameter chiptype can be used to grep files with a certain substring.
smart.install.scavenge.web.packages("http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/latest/ense.asp", chiptype = "ex", update = 1)
smart.install.scavenge.web.packages("http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/latest/ense.asp", chiptype = "hta", update = 1)

# Test if R.script has been linked to appropriate custom_cdf packages. Not perfect, but better than nothing.
# NOTE: below commands should match those listed above
db.custom.packages <- smart.install.scavenge.web.packages("http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/latest/entrezg.asp", list.only = 1)
db.custom.packages <- c(db.custom.packages, smart.install.scavenge.web.packages("http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/latest/tairg.asp", list.only = 1))
db.custom.packages <- c(db.custom.packages, smart.install.scavenge.web.packages("http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/latest/ense.asp", chiptype = "ex", list.only = 1))
db.custom.packages <- c(db.custom.packages, smart.install.scavenge.web.packages("http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/latest/ense.asp", chiptype = "hta", list.only = 1))
check.affy.customnames(script.basename, db.custom.packages)
