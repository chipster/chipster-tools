# TOOL metabarcoding-rarefy-even.R: "Rarefy OTU data to even depth" (Resamples a phyloseq OTU table such that all samples have the same sequencing depth, and saves the resulting object as an Rda file. Requires a phyloseq object in Rda format as the input.)
# INPUT ps.Rda: "Phyloseq object in Rda format" TYPE GENERIC
# OUTPUT ps_rarefied.txt: ps_rarefied.txt
# OUTPUT ps_rarefied.Rda: ps_rarefied.Rda
# RUNTIME R-3.6.1-phyloseq

# JH 2020

# Load phyloseq
library(phyloseq)

# Load phyloseq object
load("ps.Rda")

# Rarefy to even depth (after setting seed to 1)
set.seed(1)
ps <- rarefy_even_depth(ps)

# Print out object summary
sink("ps_rarefied.txt")
	cat("\n\n\n")
	cat("### Phyloseq object (after rarefying to even depth) ###\n")
	cat("\n\n\n")
	print(ps)
	cat("\n\n\n")
sink()

# Export phyloseq object as Rda file 
save(ps, file = "ps_rarefied.Rda")
