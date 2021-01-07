# TOOL metabarcoding-single-double.R: "Remove OTUs with 0-2 occurrences" (Filters a phyloseq object so that only OTUs with \>0, \>1 or \>2 occurrences are retained. The resulting object is saved in Rda format. Requires a phyloseq object in Rda format as the input.) 
# INPUT ps.Rda: "Phyloseq object in Rda format" TYPE GENERIC
# OUTPUT ps_pruned.Rda: ps_pruned.Rda
# OUTPUT ps_pruned.txt: ps_pruned.txt
# PARAMETER type: "Filtering threshold (no. of occurrences)" TYPE [0: "\>0", 1: "\>1", 2: "\>2"] DEFAULT 2 (Choose an option to retain OTUs with more than zero, single or two occurrences; the default is >2)
# RUNTIME R-3.6.1-phyloseq

# JH 2020

# Load phyloseq
library(phyloseq)

# Load phyloseq object
load("ps.Rda")

# Filter the data set
if (type == "0"){
	ps <- prune_taxa(taxa_sums(ps) > 0, ps)
	} else if (type == "1"){
	ps <- prune_taxa(taxa_sums(ps) > 1, ps)
	} else if (type == "2"){
	ps <- prune_taxa(taxa_sums(ps) > 2, ps)
}

# Print out a summary of the phyloseq object
sink("ps_pruned.txt")
	cat("\n\n\n")
	cat("### Phyloseq object ###\n")
	cat("### (following removal of OTUs with 0-2 occurrences) ###\n")
	cat("\n\n\n")
	print(ps)
	cat("\n\n\n")
sink()

# Export phyloseq object as Rda file 
save(ps, file = "ps_pruned.Rda")