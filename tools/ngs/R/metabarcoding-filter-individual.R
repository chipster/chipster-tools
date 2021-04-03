# TOOL metabarcoding-filter-individual.R: "Remove selected taxa" (This tool can be used to remove chloroplast and mitochondrial sequences from a phyloseq object, and up to five user-specified taxa at the desired level of biological organization. Produces a taxonomy summary of the resulting phyloseq object and saves it as an Rda file. Requires a phyloseq object in Rda format as the input.) 
# INPUT ps.Rda: "Phyloseq object in .Rda format" TYPE GENERIC
# OUTPUT ps_ind.Rda
# OUTPUT ps_ind_taxon.txt
# PARAMETER OPTIONAL remove.chloroplast: "Remove class Chloroplast" TYPE [yes, no] DEFAULT yes (Remove class Chloroplast)
# PARAMETER OPTIONAL remove.mitochondria: "Remove family Mitochondria" TYPE [yes, no] DEFAULT yes (Remove family Mitochondria)
# PARAMETER OPTIONAL type: "Level of biological organization for manual taxon removal" TYPE [phylum: "Phylum", class: "Class", order: "Order", family: "Family", genus: "Genus", species: "Species"] DEFAULT phylum (Select the desired taxonomic level; default is phylum)
# PARAMETER OPTIONAL tax1: "Taxon to be removed" TYPE STRING (Name of taxon to be filtered out)
# PARAMETER OPTIONAL tax2: "2nd taxon to be removed" TYPE STRING (Name of taxon to be filtered out)
# PARAMETER OPTIONAL tax3: "3rd taxon to be removed" TYPE STRING (Name of taxon to be filtered out)
# PARAMETER OPTIONAL tax4: "4th taxon to be removed" TYPE STRING (Name of taxon to be filtered out)
# PARAMETER OPTIONAL tax5: "5th taxon to be removed" TYPE STRING (Name of taxon to be filtered out)
# RUNTIME R-3.6.1-phyloseq

# JH 2020-2021

# Load packages
library(phyloseq)

# Load phyloseq object
load("ps.Rda")

# If "Species" selected for taxon removal, check that using 7-level classification
taxlength <- length(colnames(tax_table(ps)))
if (type == "species" && 
	taxlength == "6"){
		stop("Data set lacks species-level assignments; please choose another level of organization.")
}

# Remove class Chloroplast
if (remove.chloroplast == "yes") {
	ps <- subset_taxa(ps, Class != "Chloroplast")
}

# Remove family Mitochondria
<<<<<<< HEAD
if (remove.mitochondria == "yes") {
=======
if (remove.chloroplast == "yes") {
>>>>>>> master
	ps <- subset_taxa(ps, Family != "Mitochondria")
}

# Create list of taxa to filter out
# (Note: Minimap2 tool in Alignment has examples)

filterTax <- as.character(tax1)

if ( nchar(tax2) > 0 ){
	filterTax <- c(filterTax, tax2)
}
if ( nchar(tax3) > 0 ){
	filterTax <- c(filterTax, tax3)
}
if ( nchar(tax4) > 0 ){
	filterTax <- c(filterTax, tax4)
}
if ( nchar(tax5) > 0 ){
	filterTax <- c(filterTax, tax5)
}

# Filter out the specified taxa
if (type == "phylum"){
		ps <- subset_taxa(ps, !Phylum %in% filterTax)
	} else if (type == "class"){
		ps <- subset_taxa(ps, !Class %in% filterTax)
	} else if (type == "order"){
		ps <- subset_taxa(ps, !Order %in% filterTax)
	} else if (type == "family"){
		ps <- subset_taxa(ps, !Family %in% filterTax)
	} else if (type == "genus"){
		ps <- subset_taxa(ps, !Genus %in% filterTax)
	} else if (type == "species"){
		ps <- subset_taxa(ps, !Species %in% filterTax)
}

# Create taxon table following filtering
if (type == "phylum"){
		taxonsummary <- table(tax_table(ps)[, "Phylum"], exclude = NULL)
	} else if (type == "class"){
		taxonsummary <- table(tax_table(ps)[, "Class"], exclude = NULL)
	} else if (type == "order"){
		taxonsummary <- table(tax_table(ps)[, "Order"], exclude = NULL)
	} else if (type == "family"){
		taxonsummary <- table(tax_table(ps)[, "Family"], exclude = NULL)
	} else if (type == "genus"){
		taxonsummary <- table(tax_table(ps)[, "Genus"], exclude = NULL)
	} else if (type == "species"){
		taxonsummary <- table(tax_table(ps)[, "Species"], exclude = NULL)
}

# Print out taxon table
sink("ps_ind_taxon.txt")
	cat("\n\n\n")
	cat("### Summary of taxon frequencies ###\n")
	cat("### (After removing individual taxa) ###\n")
	cat("\n\n\n")
	print(taxonsummary)
	cat("\n\n\n")
sink()

# Export phyloseq object as Rda file 
save(ps, file = "ps_ind.Rda")
