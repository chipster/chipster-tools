# TOOL acgh-add-cytobands.R: "Add cytogenetic bands" (Adds the cytogenetic band information using chromosome names and start end base pair positions. If this position information is not present in your data set, please first run the Fetch probe positions from GEO/CanGEM tool.)
# INPUT normalized.tsv: normalized.tsv TYPE GENERIC 
# OUTPUT cytobands.tsv: cytobands.tsv 
# PARAMETER genome.build: "Genome build" TYPE [GRCh38: GRCh38, GRCh37: GRCh37, NCBI36: NCBI36] DEFAULT GRCh37 (The genome build to use. GRCh38 = hg38, GRCh37 = hg19, NCBI36 = hg18.)

# Ilari Scheinin <firstname.lastname@gmail.com>
# 2016-02-23

source(file.path(chipster.common.path, "library-Chipster.R"))

dat <- readData("normalized.tsv")

# check we have necessary position information
pos <- c("chromosome", "start", "end")
if (length(setdiff(pos, colnames(dat))) != 0)
  stop("CHIPSTER-NOTE: This script can only be run on files that have the following columns: chromosome, start, end.")

# change chromosome to character even if it contains only numbers
dat$chromosome <- as.character(dat$chromosome)

# load cytobands
cytodir <- list.files(file.path(chipster.tools.path, "genomes", "genomebrowser", "Homo_sapiens"), pattern=paste0("^", genome.build, "\\..*$"), full.name=TRUE)
if (length(cytodir) == 0)
  stop("CHIPSTER-NOTE: Cytoband directory not found.")
if (length(cytodir) > 1)
  cytodir <- sort(cytodir, decreasing=TRUE)[1]
cytofile <- list.files(cytodir, pattern=paste0("Homo_sapiens\\.", genome.build, "\\..*\\.cytoband-chr.txt"), full.names=TRUE)
if (length(cytofile) == 0)
  stop("CHIPSTER-NOTE: Cytoband file not found.")
if (length(cytofile) > 1)
  cytofile <- sort(cytofile, decreasing=TRUE)[1]
bands <- read.table(cytofile, sep="\t", as.is=TRUE, col.names=c("chromosome", "index", "start", "end", "band", "dye"))
bands$chromosome <- sub("^chr", "", bands$chromosome)
bands$band <- paste0(bands$chromosome, bands$band)
rownames(bands) <- bands$band
bands <- bands[bands$chromosome %in% unique(dat$chromosome), ]
bands <- bands[order(bands$index), ]

# add cytoband column after position columns
dat2 <- dat[, pos]
dat2$cytoband <- NA_character_
dat2 <- cbind(dat2, dat[, setdiff(colnames(dat), pos), drop=FALSE])
dat <- dat2
dat2 <- NULL

# add cytoband information for start and end positions
dat$startband <- "unknown"
dat$endband <- "unknown"
for (band in rownames(bands)) {
  index <- !is.na(dat$chromosome) &
    dat$chromosome == bands[band, "chromosome"] &
    dat$start      >= bands[band, "start"] &
    dat$start      <= bands[band, "end"]
  if (length(index) > 0)
    dat[index, "startband"] <- bands[band, "band"]
  index <- !is.na(dat$chromosome) &
    dat$chromosome == bands[band, "chromosome"] &
    dat$end        >= bands[band, "start"] &
    dat$end        <= bands[band, "end"]
  if (length(index) > 0)
    dat[index, "endband"] <- bands[band, "band"]
}

# format output as a range of bands or just a single band
dat$cytoband <- paste0(dat$startband, "-", dat$endband)
dat$cytoband[dat$startband == dat$endband] <- dat$startband[dat$startband == dat$endband]

# remove startband end endband from output
dat$startband <- NULL
dat$endband <- NULL

writeData(dat, "cytobands.tsv")

# EOF
