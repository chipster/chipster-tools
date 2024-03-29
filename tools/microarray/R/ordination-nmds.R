# TOOL ordination-nmds.R: NMDS (Non-metric multidimensional scaling. Creates a 2-dimensional representation of the arrays using Euclidean distances. Can be used for quality control.)
# INPUT normalized.tsv: normalized.tsv TYPE GENE_EXPRS
# INPUT META phenodata.tsv: phenodata.tsv TYPE GENERIC
# OUTPUT nmds.pdf: nmds.pdf
# PARAMETER image.width: "Image width" TYPE INTEGER FROM 200 TO 3200 DEFAULT 600 (Width of the plotted network image)
# PARAMETER image.height: "Image height" TYPE INTEGER FROM 200 TO 3200 DEFAULT 600 (Height of the plotted network image)


# Non-metric multidimensional scaling
# JTT 14.1.2008

# MG, 15.9.2011
# Updated colors and legend

# Parameter settings (default) for testing purposes
# image.width<-c(600)
# image.height<-c(600)

# Renaming variables
w <- image.width
h <- image.height

# Loads the libraries
library(MASS)

# Loads the normalized data
file <- c("normalized.tsv")
dat <- read.table(file, header = T, sep = "\t", row.names = 1)

# Separates expression values and flags
calls <- dat[, grep("flag", names(dat))]
dat2 <- dat[, grep("chip", names(dat))]

# Loads phenodata
phenodata <- read.table("phenodata.tsv", header = T, sep = "\t")

# Calculating the distance matrix
dat2.dist <- dist(t(dat2))

# Calculating the MDS
mds <- isoMDS(dat2.dist)


# Setup sample colors according to group
if (length(levels(as.factor(phenodata$group))) > 0) {
    sample_colors <- numeric(length(phenodata$group))
    group_levels <- levels(as.factor(phenodata$group))
    group_identity <- as.character(phenodata$group)
    for (count_levels in 1:length(group_levels)) {
        for (count_samples in 1:length(phenodata$group)) {
            if (group_identity[count_samples] == group_levels[count_levels]) sample_colors[count_samples] <- 1 + count_levels
        }
    }
    level_colors <- levels(as.factor(sample_colors))
}
if (length(levels(as.factor(phenodata$group))) == 0) {
    sample_colors <- rep(2, length(phenodata$group))
}


# Plotting the image
pdf(file = "nmds.pdf", width = w / 72, height = h / 72)
plot(mds$points[, 1], mds$points[, 2], main = "NMDS", pch = 19, xlab = "Dimension 1", ylab = "Dimension 2", type = "n")
text(mds$points[, 1], mds$points[, 2], phenodata$description, cex = 0.75, col = sample_colors)
if (length(levels(as.factor(phenodata$group))) > 0) {
    legend(x = "topleft", legend = group_levels, col = level_colors, cex = 0.5, pch = 19)
}
dev.off()
