# TOOL calculate-descriptive-statistics.R: "Calculate descriptive statistics" (Calculates basic descriptive statistics for all genes. These include parametric and non-parametric location and spread descriptives.)
# INPUT normalized.tsv: normalized.tsv TYPE GENE_EXPRS
# OUTPUT descr-stats.tsv: descr-stats.tsv
# OUTPUT descriptives.tsv: descriptives.tsv
# PARAMETER calculate.descriptives.for: "Calculate descriptives for" TYPE [genes: genes, chips: chips] DEFAULT genes (Descriptive statistics are calculated for...)


# Two-group parametric and non-parametric tests
# JTT 31.1.2008
# IS 1.20.2010: modified behavior for chips: now calculates the statistics for all columns

# Parameter settings (default) for testing purposes
# calculate.descriptives.for<-c("chips")

# Loads the libraries
library(genefilter)

# Loads the normalized data
file <- c("normalized.tsv")
dat <- read.table(file, header = T, sep = "\t", row.names = 1)

if (calculate.descriptives.for == "chips") {
   dat2 <- t(dat[, sapply(dat, is.numeric)])
} else {
   dat2 <- dat[, grep("^chip\\.", names(dat))]
}

# Calculates the descriptive statistics
nc <- ncol(dat2)
rsum <- rowSums(dat2)
rmean <- rsum / nc
rsd <- rowSds(dat2)
rcv <- rsd / rmean
rstexp <- rmean / rsd
rmin <- apply(X = dat2, MARGIN = 1, FUN = min)
rmax <- apply(X = dat2, MARGIN = 1, FUN = max)
rrange <- rmax - rmin
rmedian <- apply(X = dat2, MARGIN = 1, FUN = median)
riqr <- apply(X = dat2, MARGIN = 1, FUN = IQR)

# Saving the results
if (calculate.descriptives.for == "chips") {
   write.table(data.frame(average = rmean, median = rmedian, sd = rsd, cv = rcv, st.exp = rstexp, min = rmin, max = rmax, range = rrange, iqr = riqr), file = "descr-stats.tsv", sep = "\t", row.names = T, col.names = T, quote = F)
   write.table(data.frame(chip.average = rmean, chip.median = rmedian, chip.sd = rsd, chip.cv = rcv, chip.st.exp = rstexp, chip.min = rmin, chip.max = rmax, chip.range = rrange, chip.iqr = riqr), file = "descriptives.tsv", sep = "\t", row.names = T, col.names = T, quote = F)
} else {
   write.table(data.frame(dat, average = rmean, median = rmedian, sd = rsd, cv = rcv, st.exp = rstexp, min = rmin, max = rmax, range = rrange, iqr = riqr), file = "descr-stats.tsv", sep = "\t", row.names = T, col.names = T, quote = F)
   write.table(data.frame(gene.average = rmean, gene.median = rmedian, gene.sd = rsd, gene.cv = rcv, gene.st.exp = rstexp, gene.min = rmin, gene.max = rmax, gene.range = rrange, gene.iqr = riqr), file = "descriptives.tsv", sep = "\t", row.names = T, col.names = T, quote = F)
}

# EOF
