# TOOL plot-volcano-data-exists.R: "Volcano plot" (Plots a Volcano plot from results of a statistical test.)
# INPUT normalized.tsv: normalized.tsv TYPE GENE_EXPRS
# OUTPUT volcanoP.pdf: volcanoP.pdf
# PARAMETER fold.change.column: "Fold change column" TYPE COLUMN_SEL DEFAULT EMPTY (Column that contains the fold change values)
# PARAMETER p.value.column: "p-value column" TYPE COLUMN_SEL DEFAULT EMPTY (Column that contains the p-values)
# PARAMETER p.value.threshold: "p-value threshold" TYPE DECIMAL FROM 0 TO 1 DEFAULT 0.05 (P-value cut-off for significant results. Significant results are plotted with red dots)
# PARAMETER OPTIONAL gene.names: "Gene names" TYPE INTEGER FROM 0 TO 50 DEFAULT 10 (Choose to add labels to top genes on plot, select a number between 0 and 50)
# PARAMETER OPTIONAL gene.labels: "Gene labels" TYPE COLUMN_SEL DEFAULT EMPTY (Select which column for gene labels. This is not applicable if you do not select the number of genes you would like displayed on your plot.)
# PARAMETER OPTIONAL image.width: "Image width" TYPE INTEGER FROM 200 TO 3200 DEFAULT 600 (Width of the plotted network image)
# PARAMETER OPTIONAL image.height: "Image height" TYPE INTEGER FROM 200 TO 3200 DEFAULT 600 (Height of the plotted network image)

# JTT 08.11.2007: Volcano plot from existing results
# MG: 21.09.2009: Modified
# MK: 10.10.2013: Added p-value coloring

# Renaming variables
w <- image.width
h <- image.height

# Loads the normalized data
file <- c("normalized.tsv")
dat <- read.table(file, header = T, sep = "\t", row.names = 1)

# Sanity checks
if (fold.change.column == "EMPTY") {
   stop("CHIPSTER-NOTE: You haven't selected a column for fold change! Tool cannot be executed.")
}
if (p.value.column == "EMPTY") {
   stop("CHIPSTER-NOTE: You haven't selected a column for p-value! Tool cannot be executed.")
}
if (sum(grep(("p.adjusted|padj|FDR"), colnames(dat))) == 0) {
   stop("CHIPSTER-NOTE: You don't have any P-value columns in the dataset! Please run some statistical test first.")
}

# Extracts the data
expression <- dat[, grep(fold.change.column, colnames(dat))]
pvalues <- dat[, grep(p.value.column, colnames(dat))]
dot_colors <- rep(1, length(pvalues))
dot_colors[pvalues <= p.value.threshold] <- 2

# Plot label cutoff
cutoff <- sort(pvalues)[gene.names]
sign.genes <- which(pvalues <= cutoff)

labels <- row.names(dat)[sign.genes]
if (gene.labels != "EMPTY" && gene.labels != "") {
   labels <- dat[, grep(gene.labels, colnames(dat))][sign.genes]
}

# Plotting
pdf(file = "volcanoP.pdf", width = w / 72, height = h / 72)
plot(expression, -log10(pvalues), xlim = c(-max(abs(expression)), max(abs(expression))), main = "Volcano plot", pch = 19, col = dot_colors, xlab = "log2 (fold change)", ylab = "-log10 (p)", cex = 0.25)
text(x = expression[sign.genes], y = -log10(pvalues[sign.genes]), label = labels, cex = 0.5, col = "blue")
abline(h = -log10(p.value.threshold), lty = 2)
dev.off()
