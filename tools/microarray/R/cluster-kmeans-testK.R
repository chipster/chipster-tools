# TOOL cluster-kmeans-testK.R: "K-Means - estimate K" (K-means clustering of genes. Divides the genes in the selected data set into a specified number of clusters. Tests a number of different Ks and gives back a report.)
# INPUT normalized.tsv: normalized.tsv TYPE GENE_EXPRS
# OUTPUT kmeans-test.pdf: kmeans-test.pdf
# PARAMETER image.width: "Image width" TYPE INTEGER FROM 200 TO 3200 DEFAULT 600 (Width of the resampling image)
# PARAMETER image.height: "Image height" TYPE INTEGER FROM 200 TO 3200 DEFAULT 600 (Height of the resampling image)


# K-means clustering
# JTT 6.6.2008

# Parameter settings (default) for testing purposes
# image.width<-600
# image.height<-600

# Renaming variables
w <- image.width
h <- image.height

# Loads the data file
file <- c("normalized.tsv")
dat <- read.table(file, header = T, sep = "\t", row.names = 1)

# Separates expression values and flags
calls <- dat[, grep("flag", names(dat))]
dat2 <- dat[, grep("chip", names(dat))]

# Calculates the K-means clustering result with a number of different K
kmax <- c(100)
if (nrow(dat2) < 100) {
   kmax <- nrow(dat2)
}
km <- rep(NA, (kmax - 1))
i <- c(2)
while (i < kmax) {
   km[i] <- sum(kmeans(dat2, i, iter.max = 20000, nstart = 10)$withinss)
   if (i >= 3 & km[i - 1] / km[i] <= 1.01) {
      i <- kmax
   } else {
      i <- i + 1
   }
}

# Plotting the clustering
pdf(file = "kmeans-test.pdf", width = w / 72, height = h / 72)
plot(2:kmax, km, xlab = "K", ylab = "sum(withinss)", type = "b", pch = "+", main = "Terminated when change less than 1%")
dev.off()
