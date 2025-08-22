# TOOL metabarcoding-global-permanova.R: "Global PERMANOVA for OTU abundance data" (Performs a global permutational analysis of variance \(PERMANOVA, 999 iterations\) using a distance matrices derived from OTU abundance data. Accepts up to three phenodata variables. Depending on the study design, it is possible to select between a model including only main effects vs. one including main effects and all interactions. Requires an Rda file \\(ps_dist.Rda\) produced by the "Distance matrices and ordinations" tool as the input.)
# INPUT ps.Rda: "Data set in Rda format" TYPE GENERIC
# INPUT META phenodata.tsv: "Phenodata" TYPE GENERIC
# OUTPUT global_permanova.txt: global_permanova.txt
# PARAMETER OPTIONAL type: "Main effects only vs. main effects and interactions?" TYPE [main: "Main effects only", all: "Main effects and interactions"] DEFAULT main (PERMANOVA formula to be used. When interactions are included, the phenodata variables are added sequentially, and their order influences the results.)
# PARAMETER pheno1: "Phenodata variable 1" TYPE METACOLUMN_SEL (First phenodata variable used to specify the formula to be tested)
# PARAMETER OPTIONAL pheno2: "Phenodata variable 2" TYPE METACOLUMN_SEL DEFAULT EMPTY (Second phenodata variable used to specify the formula to be tested)
# PARAMETER OPTIONAL pheno3: "Phenodata variable 3" TYPE METACOLUMN_SEL DEFAULT EMPTY (Third phenodata variable used to specify the formula to be tested)
# RUNTIME R-4.4.3-phyloseq
# TOOLS_BIN ""

# JH 2020
# HJ 6.8.2025 Changed adonis -> adonis2 (adonis will probably be deprecated in the future). Added by = "margin" to main effects only tests. Main effects + interactions are by = "terms", so the order of variables matters for the main effect results.

# Load packages
library(phyloseq)
library(vegan)

# Load phyloseq object
load("ps.Rda")

# Global PERMANOVA
if (pheno2 == "EMPTY" &&
    pheno3 == "EMPTY") {
    set.seed(1)
    ps_adonis <- adonis2(ps_dist ~ get(pheno1), by = "margin", data = ps_df)
}
if (pheno2 != "EMPTY" &&
    pheno3 == "EMPTY" &&
    type == "main") {
    set.seed(1)
    ps_adonis <- adonis2(ps_dist ~ get(pheno1) + get(pheno2), by = "margin", data = ps_df)
}
if (pheno2 != "EMPTY" &&
    pheno3 == "EMPTY" &&
    type == "all") {
    set.seed(1)
    ps_adonis <- adonis2(ps_dist ~ get(pheno1) * get(pheno2), by = "terms", data = ps_df)
}
if (pheno2 != "EMPTY" &&
    pheno3 != "EMPTY" &&
    type == "main") {
    set.seed(1)
    ps_adonis <- adonis2(ps_dist ~ get(pheno1) + get(pheno2) + get(pheno3), by = "margin", data = ps_df)
}
if (pheno2 != "EMPTY" &&
    pheno3 != "EMPTY" &&
    type == "all") {
    set.seed(1)
    ps_adonis <- adonis2(ps_dist ~ get(pheno1) * get(pheno2) * get(pheno3), by = "terms", data = ps_df)
}

# Print results table
sink("global_permanova.txt")
cat("\n\n\n")
cat("### Global PERMANOVA summary ###\n")
cat("\n\n\n")
print(ps_adonis)
cat("\n\n\n")
sink()
