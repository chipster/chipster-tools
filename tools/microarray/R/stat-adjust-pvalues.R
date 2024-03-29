# TOOL stat-adjust-pvalues.R: "Adjust p-values" (Adjusts raw p-values in the selected column for multiple testing using a specified method.)
# INPUT normalized.tsv: normalized.tsv TYPE GENE_EXPRS
# INPUT META phenodata.tsv: phenodata.tsv TYPE GENERIC
# OUTPUT adjusted-p-values.tsv: adjusted-p-values.tsv
# PARAMETER column: "Column" TYPE COLUMN_SEL DEFAULT EMPTY (Data file column containing the p-values to adjust)
# PARAMETER p.value.adjustment.method: "p-value adjustment method" TYPE [Bonferroni: Bonferroni, Holm: Holm, Hochberg: Hochberg, BH: BH, BY: BY, Storey-Q: Storey-Q] DEFAULT BH (Multiple testing correction method)

# JTT 12.12.2008: P-value adjustment created

# Loads the libraries
library(multtest)
library(qvalue)

# Renaming variables
adj.method <- p.value.adjustment.method

# Loads the normalized data
file <- c("normalized.tsv")
dat <- read.table(file, header = T, sep = "\t", row.names = 1)

# Extracts raw p-values
rawp <- dat[, grep(column, colnames(dat))]

# Correction
if (adj.method == "Storey-Q") {
    adjp <- qvalue(rawp)
    adjp2 <- adjp$qvalues
} else {
    if (!(adj.method %in% c("BH", "BY"))) {
        adj.method <- tolower(adj.method)
    }
    adjp2 <- p.adjust(rawp, method = adj.method)
}

# Small manipulations
adjp2df <- as.data.frame(adjp2)
colnames(adjp2df) <- paste(colnames(dat)[grep(column, colnames(dat))], ".", adj.method, sep = "")

# Writing out data
write.table(data.frame(dat, adjp2df), file = "adjusted-p-values.tsv", sep = "\t", row.names = T, col.names = T, quote = F)
