# TOOL metabarcoding-single-double.R: "Remove OTUs with 0-20 occurrences" (Filters a phyloseq object so that only OTUs with \>0-20 occurrences are retained. The resulting object is saved in Rda format. Requires a phyloseq object in Rda format as the input.)
# INPUT ps.Rda: "Phyloseq object in Rda format" TYPE GENERIC
# OUTPUT ps_pruned.Rda: ps_pruned.Rda
# OUTPUT ps_pruned.txt: ps_pruned.txt
# PARAMETER thresh: "Filtering threshold (no. of occurrences)" TYPE INTEGER FROM 0 TO 20 DEFAULT 2 (Choose an option to retain OTUs with more that selected number of occurrences; the default is >2)
# RUNTIME R-4.4.3-phyloseq
# TOOLS_BIN ""

# JH 2020
# HJ 4.8.2025 Update for R-4.4.3 and allow giving any threshold from 0 to 20, instead of 0/1/2

# Load phyloseq
library(phyloseq)

# Load phyloseq object
load("ps.Rda")

# Filter the data set
ps <- prune_taxa(taxa_sums(ps) > thresh, ps)


# Print out a summary of the phyloseq object
sink("ps_pruned.txt")
cat("\n\n\n")
cat("### Phyloseq object ###\n")
cat("### (following removal of OTUs with 0-20 occurrences) ###\n")
cat("\n\n\n")
print(ps)
cat("\n\n\n")
sink()

# Export phyloseq object as Rda file
save(ps, file = "ps_pruned.Rda")
