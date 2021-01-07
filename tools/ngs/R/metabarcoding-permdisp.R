# TOOL metabarcoding-permdisp.R: "PERMDISP for OTU abundance data" (Performs a permutation-based test of the multivariate homogeneity of group dispersions \(PERMDISP, 999 iterations\), using a distance matrix derived from OTU abundance data. Accepts a single phenodata variable at a time, with the option to run three consecutive tests \(for three different variables\). Requires an Rda file \(ps_dist.Rda\) produced by the "Distance matrices and ordinations" tool as the input.)
# INPUT ps.Rda: "Data set in .Rda format" TYPE GENERIC
# INPUT META phenodata.tsv: "Phenodata" TYPE GENERIC
# OUTPUT permdisp.txt: permdisp.txt
# OUTPUT ps_disp.Rda: ps_disp.Rda
# PARAMETER pheno1: "Phenodata variable 1" TYPE METACOLUMN_SEL (Phenodata variable used for first PERMDISP analysis)
# PARAMETER OPTIONAL pheno2: "Phenodata variable 2" TYPE METACOLUMN_SEL DEFAULT EMPTY (Phenodata variable used for second PERMDISP analysis)
# PARAMETER OPTIONAL pheno3: "Phenodata variable 3" TYPE METACOLUMN_SEL DEFAULT EMPTY (Phenodata variable used for third PERMDISP analysis)
# RUNTIME R-3.6.1-phyloseq

# JH 2020

# Load packages
library(phyloseq)
library(vegan)

# Load phyloseq object
load("ps.Rda")

# PERMDISP (betadisper + permutest functions)
set.seed(1)
beta1 <- betadisper(ps_dist, ps_df[, pheno1])
set.seed(1)
permutest_beta1 <- permutest(beta1)

if (pheno2 == "EMPTY" &&
	pheno3 != "EMPTY"){
		stop("3rd phenodata variable specified without specifying 2nd; if running two analyses, specify phenodata variables 1 and 2.")
}

if (pheno2 != "EMPTY"){
	set.seed(1)
	beta2 <- betadisper(ps_dist, ps_df[, pheno2])
	set.seed(1)
	permutest_beta2 <- permutest(beta2)
}
if (pheno3 != "EMPTY"){
	set.seed(1)
	beta3 <- betadisper(ps_dist, ps_df[, pheno3])
	set.seed(1)
	permutest_beta3 <- permutest(beta3)
}

# Print results table
if (pheno2 == "EMPTY" &&
	pheno3 == "EMPTY"){
sink("permdisp.txt")
	cat("\n\n\n")
	cat("### PERMDISP summary ###\n")
	cat("\n\n\n")
	print(permutest_beta1)
	cat("\n\n\n")
sink()
}
if (pheno2 != "EMPTY" &&
	pheno3 == "EMPTY"){
sink("permdisp.txt")
	cat("\n\n\n")
	cat("### PERMDISP summary ###\n")
	cat("\n\n\n")
	cat("### Test 1 ###\n")
	cat("\n\n\n")
	print(permutest_beta1)
	cat("\n\n\n")
	cat("### Test 2 ###\n")
	cat("\n\n\n")
	print(permutest_beta2)
	cat("\n\n\n")
sink()
}
if (pheno2 != "EMPTY" &&
	pheno3 != "EMPTY"){
sink("permdisp.txt")
	cat("\n\n\n")
	cat("### PERMDISP summary ###\n")
	cat("\n\n\n")
	cat("### Test 1 ###\n")
	cat("\n\n\n")
	print(permutest_beta1)
	cat("\n\n\n")
	cat("### Test 2 ###\n")
	cat("\n\n\n")
	print(permutest_beta2)
	cat("\n\n\n")
	cat("### Test 3 ###\n")
	cat("\n\n\n")
	print(permutest_beta3)
	cat("\n\n\n")
sink()
}

# Save the dispersion and test data as Rda
if (pheno2 == "EMPTY" &&
	pheno3 == "EMPTY"){
		save(list=c("ps_df",
			"pheno1", "beta1", "permutest_beta1"), 
				file="ps_disp.Rda")
}
if (pheno2 != "EMPTY"){
	save(list=c("ps_df", 
			"pheno1", "beta1", "permutest_beta1", 
			"pheno2", "beta2", "permutest_beta2"), 
				file="ps_disp.Rda")
}
if (pheno3 != "EMPTY"){
	save(list=c("ps_df", 
			"pheno1", "beta1", "permutest_beta1",
			"pheno2", "beta2", "permutest_beta2",
			"pheno3", "beta3", "permutest_beta3"), 
				file="ps_disp.Rda")
}
