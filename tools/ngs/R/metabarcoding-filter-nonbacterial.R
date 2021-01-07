# TOOL metabarcoding-filter-nonbacterial.R: "Filter out non-bacterial features" (Tidies a phyloseq object so that only bacterial OTUs are retained. Filters out OTUs that have not been classified as Bacteria \(kingdom level\), that have been classified as NA at the phylum level, and\/or have been classified as chloroplast sequences at the class level. Produces a phylum-level taxonomy summary and prevalence table following filtering, and saves the resulting phyloseq object as an Rda file. Requires a phyloseq object in Rda format as the input.) 
# INPUT ps.Rda: "Phyloseq object in Rda format" TYPE GENERIC
# OUTPUT ps_bacteria.Rda: ps_bacteria.Rda
# OUTPUT ps_bacteria_taxon.txt: ps_bacteria_taxon.txt 
# RUNTIME R-3.6.1-phyloseq

# JH 2020

# Load packages
library(phyloseq)
library(dplyr)
library(plyr)

# Load phyloseq object
load("ps.Rda")

# Filter out anything not classified as Bacteria at the kingdom level
# (while also removing chloroplast sequences at the Class level) 
ps <- subset_taxa(ps, Kingdom == "Bacteria" & Class != "Chloroplast")

# Filter out any remaining features with the phylum annotated as NA
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

# Export phyloseq object as Rda file 
save(ps, file = "ps_bacteria.Rda")
