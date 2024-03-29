# TOOL average-replicates.R: "Average replicate chips" (Calculates averages for replicate chips.)
# INPUT normalized.tsv: normalized.tsv TYPE GENE_EXPRS
# INPUT META phenodata.tsv: phenodata.tsv TYPE GENERIC
# OUTPUT average-replicates.tsv: average-replicates.tsv
# OUTPUT META phenodata-average.tsv: phenodata.tsv
# PARAMETER column: Column TYPE METACOLUMN_SEL DEFAULT group (Phenodata column describing the groups to average.)
# PARAMETER averaging: Averaging TYPE [mean: mean, median: median] DEFAULT mean (Averaging using the mean or median.)

# Average replicate chips
# JTT 30.7.2007
#
# modified by MG, 12.4.2010
# rewritten by IS, 6.9.2010
# modified by IS, 16.2.2011
# modified by IS,12.10.2012

# load inputs
file <- "normalized.tsv"
dat <- read.table(file, header = TRUE, sep = "\t", quote = "", row.names = 1, check.names = FALSE)
phenodata <- read.table("phenodata.tsv", header = TRUE, sep = "\t", as.is = TRUE)

# identify replicates
replicates <- unique(phenodata[duplicated(phenodata[, column]), column])
unique.chips <- which(!phenodata[, column] %in% replicates)

# generate new phenodata
concatenate.if.not.equal <- function(x) {
    x <- unique(x)
    paste(x, collapse = ";")
}
numeric.cols <- sapply(phenodata, is.numeric)
other.cols <- !sapply(phenodata, is.numeric)
phenodata2 <- phenodata[unique.chips, ]
for (s in replicates) {
    ss <- phenodata[phenodata[, column] == s, ]
    ss[1, numeric.cols] <- apply(as.data.frame(ss[, numeric.cols]), 2, averaging)
    ss[1, other.cols] <- apply(as.data.frame(ss[, other.cols]), 2, concatenate.if.not.equal)
    ss <- ss[1, ]
    phenodata2 <- rbind(phenodata2, ss)
}
phenodata2$sample <- sprintf("microarray%.3i", 1:nrow(phenodata2))

# identify matrices (chip, flag, segmented, ...) present in the data
x <- colnames(dat)
suffix <- sub("^chip\\.", "", x[grep("^chip\\.", x)[1]])
matrices <- sub(suffix, "", x[grep(suffix, x)])

# identify annotation columns (that are not part of any of the matrices)
annotations <- 1:ncol(dat)
for (m in matrices) {
    annotations <- setdiff(annotations, grep(m, x))
}
dat2 <- dat[, annotations]

# generate new data table
for (m in matrices) {
    m2 <- dat[, grep(m, x)]
    m3 <- as.data.frame(m2[, unique.chips])
    num <- is.numeric(as.matrix(m2))
    for (s in replicates) {
        ss <- which(phenodata[, column] == s)
        if (num) {
            ss <- apply(m2[, ss], 1, averaging, na.rm = TRUE) # what to do with the log transformation?
        } else {
            ss <- apply(m2[, ss], 1, concatenate.if.not.equal)
        }
        m3 <- cbind(m3, ss, stringsAsFactors = FALSE)
    }
    colnames(m3) <- paste(m, phenodata2$sample, sep = "")
    dat2 <- cbind(dat2, m3)
}

# write output
options(scipen = 10)
write.table(dat2, file = "average-replicates.tsv", quote = FALSE, sep = "\t")
write.table(phenodata2, file = "phenodata-average.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

# EOF
