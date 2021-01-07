# TOOL metabarcoding-filter-individual.R: "Filter out selected phyla" (Removes up to five user-specified phyla from a phyloseq object. Produces a taxonomy summary at the phylum level and saves the resulting phyloseq object as an Rda file. Requires a phyloseq object in Rda format as the input.) 
# INPUT ps.Rda: "Phyloseq object in .Rda format" TYPE GENERIC
# OUTPUT ps_ind.Rda
# OUTPUT ps_ind_taxon.txt
# PARAMETER phylum1: "Phylum to be removed" TYPE STRING (Name of phylum to be filtered out)
# PARAMETER OPTIONAL phylum2: "2nd phylum to be removed" TYPE STRING (Name of phylum to be filtered out)
# PARAMETER OPTIONAL phylum3: "3rd phylum to be removed" TYPE STRING (Name of phylum to be filtered out)
# PARAMETER OPTIONAL phylum4: "4th phylum to be removed" TYPE STRING (Name of phylum to be filtered out)
# PARAMETER OPTIONAL phylum5: "5th phylum to be removed" TYPE STRING (Name of phylum to be filtered out)
# RUNTIME R-3.6.1-phyloseq

# JH 2020

# Load packages
library(phyloseq)

# Load phyloseq object
load("ps.Rda")

# Create list of phyla to filter out
# (Note: Minimap2 tool in Alignment has examples)

filterPhyla <- as.character(phylum1)

if ( nchar(phylum2) > 0 ){
	filterPhyla <- c(filterPhyla, phylum2)
}
if ( nchar(phylum3) > 0 ){
	filterPhyla <- c(filterPhyla, phylum3)
}
if ( nchar(phylum4) > 0 ){
	filterPhyla <- c(filterPhyla, phylum4)
}
if ( nchar(phylum5) > 0 ){
	filterPhyla <- c(filterPhyla, phylum5)
}

# Filter out the specified phyla
ps <- subset_taxa(ps, !Phylum %in% filterPhyla)

# Create phylum-level taxon table
taxonsummary <- table(tax_table(ps)[, "Phylum"], exclude = NULL)

# Print out phylum-level taxon table
sink("ps_ind_taxon.txt")
	cat("\n\n\n")
	cat("### Summary of taxon frequencies ###\n")
	cat("### (After removing invidivual phyla) ###\n")
	cat("\n\n\n")
	print(taxonsummary)
	cat("\n\n\n")
sink()

# Export phyloseq object as Rda file 
save(ps, file = "ps_ind.Rda")
