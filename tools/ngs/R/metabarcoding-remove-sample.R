# TOOL metabarcoding-remove-sample.R: "Remove samples from the Phyloseq object" (Remove one to five samples from a phyloseq object. Sample names must be given exactly as they are written in the phenodata table.)
# INPUT ps.Rda: "Phyloseq object in Rda format" TYPE GENERIC
# INPUT META phenodata.tsv: "Phenodata" TYPE GENERIC
# OUTPUT ps_subset.Rda: ps_subset.Rda
# OUTPUT ps_subset.txt: ps_subset.txt
# PARAMETER xvar: "Phenodata variable with sequencing sample IDs" TYPE METACOLUMN_SEL DEFAULT EMPTY (Phenodata variable with unique IDs for each community profile.)
# PARAMETER samplex1: "1st sample to be removed" TYPE STRING (Type the name of the sample as it appears in the phenodata)
# PARAMETER OPTIONAL samplex2: "2nd sample to be removed" TYPE STRING (Type the name of the sample as it appears in the phenodata)
# PARAMETER OPTIONAL samplex3: "3rd sample to be removed" TYPE STRING (Type the name of the sample as it appears in the phenodata)
# PARAMETER OPTIONAL samplex4: "4th sample to be removed" TYPE STRING (Type the name of the sample as it appears in the phenodata)
# PARAMETER OPTIONAL samplex5: "5th sample to be removed" TYPE STRING (Type the name of the sample as it appears in the phenodata)
# RUNTIME R-4.4.3-phyloseq
# TOOLS_BIN ""

# HJ 2025

# Load packages
library(phyloseq)
library(ggplot2)

# Load phyloseq object
load("ps.Rda")

# create a vector of samples to be removed
removeSamples <- as.character(samplex1)

if (nchar(samplex2) > 0) {
    removeSamples <- c(removeSamples, samplex2)
}
if (nchar(samplex3) > 0) {
    removeSamples <- c(removeSamples, samplex3)
}
if (nchar(samplex4) > 0) {
    removeSamples <- c(removeSamples, samplex4)
}
if (nchar(samplex5) > 0) {
    removeSamples <- c(removeSamples, samplex5)
}

# check that the samples to be removed exist in the phyloseq object and give warning if not
#if (removeSamples %in% sample_names(ps) == FALSE) {
if (length(intersect(removeSamples, sample_names(ps))) != length(removeSamples)) {
        stop("CHIPSTER-NOTE: A sample to be removed is not found in the Phyloseq object. Please check the sample names.")
    }

# remove the given samples from the phyloseq object
#ps = subset_samples(ps, sample != removeSamples)
ps = subset_samples(ps, !sample %in% removeSamples)

sink("ps_subset.txt")
cat("\n\n\n")
cat("### Summary of Phyloseq object ###\n")
cat("### (After removal of specified samples) ###\n")
cat("\n\n\n")
print(ps)
cat("\n\n\n")
sink()

# save the subsetted phyloseq object
save(ps, file = "ps_subset.Rda")
