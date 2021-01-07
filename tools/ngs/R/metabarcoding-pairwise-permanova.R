# TOOL metabarcoding-pairwise-permanova.R: "Post-hoc pairwise PERMANOVA for OTU abundance data" (Performs post-hoc pairwise permutational analyses of variance \(PERMANOVA, 999 iterations\), using distance matrices derived from OTU abundance data. Accepts a single phenodata variable and applies a Benjamini-Hochberg correction to account for multiple testing. Requires an Rda file \(ps_dist.Rda\) produced by the "Distance matrices and ordinations" tool as the input. This tool should only be used after obtaining a statistically significant \(p \< 0.05\) result using the \"Global PERMANOVA for OTU abundance data\" tool.)
# INPUT ps.Rda: "Data set in Rda format" TYPE GENERIC
# INPUT META phenodata.tsv: "Phenodata" TYPE GENERIC
# OUTPUT pairwise_permanova.txt: pairwise_permanova.txt
# PARAMETER pheno: "Phenodata variable" TYPE METACOLUMN_SEL (Phenodata variable used for pairwise comparisons)
# RUNTIME R-3.6.1-phyloseq

# JH 2020

# Load packages
library(phyloseq)
library(RVAideMemoire)

# Load phyloseq object
load("ps.Rda")

# Pairwise PERMANOVA
set.seed(1)
ps_pairwise <- pairwise.perm.manova(ps_dist, ps_df[, pheno], 
	nperm = 999, p.method = "BH")

# Print results table
sink("pairwise_permanova.txt")
	cat("\n\n\n")
	cat("### Pairwise PERMANOVA summary ###\n")
	cat("\n\n\n")
	print(ps_pairwise)
	cat("\n\n\n")
sink()
