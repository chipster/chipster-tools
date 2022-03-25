# TOOL ordination-pca.R: PCA (Principal component analysis. The number of principal component to save is controlled through the explained variablity. All principal components are saved until the explained variability is exceeded, but at least 3 components are always saved.)
# INPUT normalized.tsv: normalized.tsv TYPE GENE_EXPRS 
# INPUT META phenodata-normalized.tsv: "Phenodata" TYPE GENERIC
# OUTPUT pca.tsv: pca.tsv 
# OUTPUT variance.tsv: variance.tsv 
# OUTPUT loadings.tsv: loadings.tsv
# OUTPUT pca_biplot.pdf: pca_biplot.pdf
# PARAMETER samplevar: "Phenodata variable with sequencing sample IDs" TYPE METACOLUMN_SEL DEFAULT EMPTY (Phenodata variable with unique ID for each sample)
# PARAMETER OPTIONAL group_column: "Phenodata variable for PCA grouping" TYPE METACOLUMN_SEL DEFAULT empty (Phenodata variable describing grouping added to PCA biplot)
# PARAMETER OPTIONAL expvar: "Amount of variation to explain" TYPE INTEGER FROM 0 TO 100 DEFAULT 80 (Percentage of experimental variation to explain, as a percentage 0-100)
# PARAMETER OPTIONAL scaling: Scaling TYPE [yes: yes, no: no] DEFAULT no (Scale the data to have a unit variance)
# PARAMETER OPTIONAL centering: Centering TYPE [yes: yes, no: no] DEFAULT yes (Scale the data to have the same mean)
# PARAMETER OPTIONAL vectors: "No. of top genes to plot as vectors" TYPE INTEGER FROM 0 TO 5 DEFAULT 0 (Number of top genes to plot as vectors in PCA \(max 5\))
# IMAGE comp-20.04-r-deps
# RUNTIME R-4.1.1-statistics

# JTT, 27.6.2006
# MG, 29.3.2010 to enable coloring by cluster feature 
# IS, 12.10.2012 to cope with tables with gene descriptions (that typically contain 's)
# MK, 16.09.2013 modified to produce variance explanation table
# EK, 02.07.2014 modified to produce loadings table
# JH, 2022 modified for PCA biplot functionality, construction of normalized + phenodata joint table

# Load factoextra git
library(factoextra)

# Load normalized data
file <- 'normalized.tsv'
dat <- read.table(file, header = TRUE,
                  sep = '\t', quote = '', row.names = 1, 
                  check.names=FALSE)

# Load phenodata
phenodata <- read.delim("phenodata-normalized.tsv")

# Separate expression values and other columns
dat2 <- dat[,grep("chip", names(dat))]

# Transpose the data
dat2 <- t(dat2)

# Remove "chip." from beginning of dat2 rownames
# (Otherwise Chipster automatically produces scatter + expression plots that
# are not required for this particular tool)

rownames(dat2) <- gsub("chip.", "", rownames(dat2))

# PCA calculations
if(scaling == "yes" & centering == "yes") {
  pc <- prcomp(dat2, scale = T, center = T)
}
if(scaling=="yes" & centering=="no") {
  pc <- prcomp(dat2, scale = T, center = F)
}
if(scaling == "no" & centering == "yes") {
  pc <- prcomp(dat2, scale = F, center = T)
}
if(scaling == "no" & centering == "no") {
  pc <- prcomp(dat2, scale = F, center = F)
}

# How many PCs to save? If too low, three will be printed at minimum
no<-as.vector(head(which(summary(pc)$importance[3,]>=(expvar/100)), n=1)[1])
if(no<3) {
  no<-c(3)
}

# Convert PCs from matrix format into data frame
pcs<-as.data.frame(pc$x[,1:no])

# Give the PC headers new names
# MAR 22: REMOVED
# names(pcs)<-paste("chip.", names(pcs), sep="")

# Round PCs to two digits
pcs<-round(pcs,digits=2)

# Save the PCs with data
write.table(data.frame(pcs), "pca.tsv", sep="\t", row.names=T, col.names=T, quote=F)

# Variance explained
write.table(round(data.frame(summary(pc)$importance[,1:no]),3), "variance.tsv", sep="\t", row.names=T, col.names=T, quote=F)

# Output component loadings. Add ".chip" to column names so that min and max can be easily found with the tool "Calculate descriptive statistics".
loadings<-round(data.frame(pc$rotation[,1:no]),6)
colnames(loadings)<-paste("chip.", colnames(loadings), sep="")
write.table(loadings, "loadings.tsv", sep="\t", row.names=T, col.names=T, quote=F)

# Biplot preparatory steps

# Set phenodata levels using names in sample column
# to test: 
# samplevar <- "sample"

# Note: if wanted to add "chip." here, could do as follows:
# levels(phenodata[, samplevar[1]]) <- paste("chip.", phenodata[, samplevar[1]], sep = "")

levels(phenodata[, samplevar[1]]) <- paste(phenodata[, samplevar[1]], sep = "")

# Rename phenodata rownames using levels
rownames(phenodata) <- levels(phenodata[, samplevar[1]])

# Convert dat2 to data frame
dat2 <- as.data.frame(dat2)

# Merge dat2 and phenodata, remove "Row.names" column produced by merge()

# Note: why do we need dat3?
# prcomp() wants only numeric columns, specifying these in practice
# is difficult given that phenodata variables in Chipster will vary between analyses.
# A simpler approach is to use the same df for plotting, but with phenodata
# variables appended.

dat3 <- transform(merge(dat2, phenodata, by = 0, all = TRUE), 
                  row.names = Row.names, Row.names = NULL)

# Biplot construction

# to test:
# group_column <- "group"
# vectors <- "4"

vectors <- as.numeric(vectors)

# plot without vectors
if(vectors == 0) {
  pca_biplot <- fviz_pca_ind(pc, 
                             label = "none",
                             repel = TRUE,
                             habillage = dat3[, group_column[1]],
                             addEllipses = TRUE, 
                             ellipse.level = 0.95) +
    theme_minimal() +
    ggtitle("PCA")
}

# plot with vectors
# show top genes explaining differences

if(vectors > 0) {
  pca_biplot <- fviz_pca_biplot(pc,
                                label ="var",
                                repel = TRUE,
                                select.var = list(contrib = vectors),
                                habillage = dat3[, group_column[1]],
                                addEllipses = TRUE,
                                ellipse.type = "confidence",
                                ellipse.level = 0.95) +
    theme_minimal() +
    ggtitle("PCA (top contributing genes shown as vectors)")
}

# Open a report PDF
pdf("pca_biplot.pdf")

pca_biplot

# Close the report PDF
dev.off()

# EOF