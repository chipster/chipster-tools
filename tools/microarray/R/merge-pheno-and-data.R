# TOOL merge-pheno-and-data.R: "Merge expression and phenodata" (Merge phenodata and expression matrices into a single spreadsheet for easy data exporting.)
# INPUT normalized.tsv: "Expression data" TYPE GENE_EXPRS
# INPUT META phenodata.tsv: phenodata.tsv TYPE GENERIC
# OUTPUT merged.tsv
# PARAMETER excel.file: "Produce Excel file" TYPE [yes, no] DEFAULT no (Fix column headers so that the file opens correctly in Excel)

# MK 20.09.2013
# EK 19.11.2013 Text changes.
# ML 17.07.2015 Fixes to the parameter names and explanations

# Sanity checks
file <- c("phenodata.tsv")
dat1 <- read.table(file, sep = "\t", header = T, row.names = 1, quote = "")

file <- c("normalized.tsv")
dat2 <- read.table(file, sep = "\t", header = T, row.names = 1, quote = "", check.names = FALSE)

# create an empty matrix
temp <- data.frame(matrix(NA, ncol = ncol(dat2), nrow = ncol(dat1) + 1))
colnames(temp) <- colnames(dat2)
dat3 <- rbind(temp, dat2)

# find expression columns associated with the first sample and check what it their prefix
firstsample <- colnames(dat2)[grep(rownames(dat1)[1], colnames(dat2))]
prefixes <- gsub(rownames(dat1)[1], "", firstsample)

for (i in 1:length(prefixes)) {
    temp <- dat1[!is.na(match(paste(prefixes[i], rownames(dat1), sep = ""), colnames(dat2))), ]
    cols <- match(paste(prefixes[i], rownames(temp), sep = ""), colnames(dat2))
    dat3[1, cols] <- rownames(temp)
    dat3[2:(ncol(dat1) + 1), cols] <- t(temp)
}

# rownames(dat3) <- make.names(c("sample", colnames(dat1), rownames(dat2)), unique=TRUE)
dat3 <- data.frame(dat3, row.names = c("sample", colnames(dat1), rownames(dat2)), check.names = FALSE)

# Writes out the combined table
if (excel.file == "no") {
    write.table(dat3, "merged.tsv", sep = "\t", row.names = T, col.names = T, quote = F)
} else {
    write.table(dat3, "merged.tsv", sep = "\t", row.names = T, col.names = NA, quote = F)
}
