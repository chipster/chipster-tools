# TOOL metabarcoding-taxon-overview.R: "Overview of taxon composition" (Tabulates the overall OTU count and a summary of taxon frequencies. Requires a phyloseq object in Rda format as the input.)
# INPUT ps.Rda: "Phyloseq object in Rda format" TYPE GENERIC
# OUTPUT ps_taxon.txt: ps_taxon.txt
# PARAMETER OPTIONAL type: "Level of biological organization for tabulating taxon composition" TYPE [phylum: "Phylum", class: "Class", order: "Order", family: "Family", genus: "Genus", species: "Species"] DEFAULT phylum (Select the desired taxonomic level; default is phylum)
# RUNTIME R-3.6.1-phyloseq

# JH 2020

# Load packages
library(phyloseq)

# Load phyloseq object
load("ps.Rda")

# Overall no. of taxa
taxano <- ntaxa(ps)

# if "Species" selected for taxon tabulation, check that using 7-level classification
taxlength <- length(colnames(tax_table(ps)))
if (type == "species" && 
	taxlength == "6"){
		stop("Data set lacks species-level assignments; please choose another level of organization.")
}

# Taxonomy table
if (type == "phylum"){
	taxonsummary <- table(tax_table(ps)[, "Phylum"], exclude = NULL) # column no. equivalents would be e.g. [, 2]
	} else if (type == "class"){
		taxonsummary <- table(tax_table(ps)[, "Class"], exclude = NULL)
	} else if (type == "order"){
		taxonsummary <- table(tax_table(ps)[, "Order"], exclude = NULL)
	} else if (type == "family"){
		taxonsummary <- table(tax_table(ps)[, "Family"], exclude = NULL)
	} else if (type == "genus") {
		taxonsummary <- table(tax_table(ps)[, "Genus"], exclude = NULL)
	} else if (type == "species") {
		taxonsummary <- table(tax_table(ps)[, "Species"], exclude = NULL)		
}

sink("ps_taxon.txt")
	cat("\n\n\n")
	cat("### Overall no. of OTUs ###\n")
	cat("\n\n\n")
	print(taxano)
	cat("\n\n\n")
	cat("### Summary of taxon frequencies ###\n")
	cat("\n\n\n")
	print(taxonsummary)
	cat("\n\n\n")
sink()
