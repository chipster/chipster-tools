# TOOL metabarcoding-filter-taxgroup.R: "Filter by taxonomic group" (Tidies a phyloseq object so that OTUs only from the desired taxonomic group \(bacteria, archaea, eukaryotes or fungi\) are retained. If bacteria are selected as the group, filters out OTUs not classified as Bacteria \(domain level\), classified as NA at the phylum level, and\/or classified as chloroplast sequences at the class level. With the exception of chloroplast sequence removal, similar filtering steps are applied to other groups. For fungi, filtering is performed at the kingdom level \(Fungi\) rather than the domain level. Produces a phylum-level taxonomy summary and prevalence table following filtering, and saves the resulting phyloseq object as an Rda file. Requires a phyloseq object in Rda format as the input.)
# INPUT ps.Rda: "Phyloseq object in Rda format" TYPE GENERIC
# OUTPUT OPTIONAL ps_bacteria.Rda
# OUTPUT OPTIONAL ps_bacteria_taxon.txt
# OUTPUT OPTIONAL ps_archaea.Rda
# OUTPUT OPTIONAL ps_archaea_taxon.txt
# OUTPUT OPTIONAL ps_eukaryota.Rda
# OUTPUT OPTIONAL ps_eukaryota_taxon.txt
# OUTPUT OPTIONAL ps_fungi.Rda
# OUTPUT OPTIONAL ps_fungi_taxon.txt
# PARAMETER group: "Group to retain" TYPE [bacteria: "Bacteria", archaea: "Archaea", eukaryotes: "Eukaryotes", fungi: "Fungi"] DEFAULT bacteria (Taxonomic group to retain)
# RUNTIME R-3.6.1-phyloseq

# JH 2020-2021

# Load packages
library(phyloseq)
library(dplyr)
library(plyr)

# Load phyloseq object
load("ps.Rda")

# Group-specific filtering

if (group == "bacteria"){
# Filter out anything not classified as Bacteria at the domain level
# (while also removing chloroplast sequences at the Class level) 
ps <- subset_taxa(ps, Domain_Kingdom == "Bacteria" & Class != "Chloroplast")
}

if (group == "archaea"){
# Filter out anything not classified as Archaea at the domain level 
ps <- subset_taxa(ps, Domain_Kingdom == "Archaea")
}

if (group == "eukaryotes"){
# Filter out anything not classified as Eukaryota at the domain level
ps <- subset_taxa(ps, Domain_Kingdom == "Eukaryota")
}

if (group == "fungi"){
# Filter out anything not classified as Fungi at the kingdom level
ps <- subset_taxa(ps, Domain_Kingdom == "Fungi")
}

# General steps

# Remove any remaining features with the phylum annotated as NA
# (or with otherwise ambiguous phylum-level annotation)
ps <- subset_taxa(ps, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))

# Phylum-level taxonomy table
taxonsummary <- table(tax_table(ps)[, "Phylum"], exclude = NULL)

# Prevalences of each feature (stored as data frame)
prevdf <- apply(X = otu_table(ps),
               MARGIN = ifelse(taxa_are_rows(ps), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})

# Add taxonomy and total read counts to data frame
prevdf <- data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(ps),
                    tax_table(ps))

# Tabulate average and total prevalences
prevsummary <- ddply(prevdf, "Phylum", 
		function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})

# Give the columns better titles
prevsummary <- dplyr::rename(prevsummary, 
		Mean.prevalence = 2, Total.prevalence = 3) # Rename cols. 2 and 3

# Print out results

if (group == "bacteria"){
sink("ps_bacteria_taxon.txt")
	cat("\n\n\n")
	cat("### Summary of phylum frequencies ###\n")
	cat("### (Following removal of non-bacterial / chloroplast features) ###\n")
	cat("\n\n\n")
	print(taxonsummary)
	cat("\n\n\n")
	cat("### Phylum-level mean and total prevalences ###\n")
	cat("\n\n\n")
	print(prevsummary)
	cat("\n\n\n")
sink()
}

if (group == "archaea"){
	sink("ps_archaea_taxon.txt")
	cat("\n\n\n")
	cat("### Summary of phylum frequencies ###\n")
	cat("### (Following removal of non-archaeal features) ###\n")
	cat("\n\n\n")
	print(taxonsummary)
	cat("\n\n\n")
	cat("### Phylum-level mean and total prevalences ###\n")
	cat("\n\n\n")
	print(prevsummary)
	cat("\n\n\n")
sink()
}

if (group == "eukaryotes"){
	sink("ps_eukaryota_taxon.txt")
	cat("\n\n\n")
	cat("### Summary of phylum frequencies ###\n")
	cat("### (Following removal of non-eukaryotic features) ###\n")
	cat("\n\n\n")
	print(taxonsummary)
	cat("\n\n\n")
	cat("### Phylum-level mean and total prevalences ###\n")
	cat("\n\n\n")
	print(prevsummary)
	cat("\n\n\n")
sink()
}

if (group == "fungi"){
	sink("ps_fungi_taxon.txt")
	cat("\n\n\n")
	cat("### Summary of phylum frequencies ###\n")
	cat("### (Following removal of non-fungal features) ###\n")
	cat("\n\n\n")
	print(taxonsummary)
	cat("\n\n\n")
	cat("### Phylum-level mean and total prevalences ###\n")
	cat("\n\n\n")
	print(prevsummary)
	cat("\n\n\n")
sink()
}

# Export phyloseq object as Rda file 

if (group == "bacteria"){
	save(ps, file = "ps_bacteria.Rda")
}

if (group == "archaea"){
	save(ps, file = "ps_archaea.Rda")
}

if (group == "eukaryotes"){
	save(ps, file = "ps_eukaryota.Rda")
}

if (group == "fungi"){
	save(ps, file = "ps_fungi.Rda")
}
