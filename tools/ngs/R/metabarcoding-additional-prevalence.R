# TOOL metabarcoding-additional-prevalence.R: "Additional prevalence summaries" (Prevalence summaries to assist with data preparation for downstream analyses. Plots total read counts against prevalence values for each phylum. Also tabulates the overall numbers of features with zero reads, one read \(singletons\) and two reads \(doubletons\). Requires a phyloseq object in Rda format as the input.) 
# INPUT ps.Rda: "Phyloseq object in .Rda format" TYPE GENERIC
# OUTPUT ps_prevalence.pdf: ps_prevalence.pdf
# OUTPUT ps_low.txt: ps_low.txt
# RUNTIME R-3.6.1-phyloseq

# JH 2020

# Load packages
library(phyloseq)
library(data.table)
library(ggplot2)

# Load phyloseq object
load("ps.Rda")

# Define fast_melt() function
# This function is originally described in taxa_summary.R on the evomics website:
# http://evomics.org/phyloseq/

fast_melt = function(physeq,
                     includeSampleVars = character(),
                     omitZero = FALSE){
  otutab = as(otu_table(physeq), "matrix")
  if(!taxa_are_rows(physeq)){otutab <- t(otutab)}
  otudt = data.table(otutab, keep.rownames = TRUE)
  setnames(otudt, "rn", "TaxaID")
  otudt[, TaxaIDchar := as.character(TaxaID)]
  otudt[, TaxaID := NULL]
  setnames(otudt, "TaxaIDchar", "TaxaID")
  mdt = melt.data.table(otudt, 
                        id.vars = "TaxaID",
                        variable.name = "SampleID",
                        value.name = "count")
  if(omitZero){
    mdt <- mdt[count > 0]
  }
  mdt <- mdt[!is.na(count)]
  mdt[, RelativeAbundance := count / sum(count), by = SampleID]
  if(!is.null(tax_table(physeq, errorIfNULL = FALSE))){
    taxdt = data.table(as(tax_table(physeq, errorIfNULL = TRUE), "matrix"), keep.rownames = TRUE)
    setnames(taxdt, "rn", "TaxaID")
    taxdt[, TaxaIDchar := as.character(TaxaID)]
    taxdt[, TaxaID := NULL]
    setnames(taxdt, "TaxaIDchar", "TaxaID")
    setkey(taxdt, "TaxaID")
    setkey(mdt, "TaxaID")
    mdt <- taxdt[mdt]
  }
  wh.svars = which(sample_variables(physeq) %in% includeSampleVars)
  if( length(wh.svars) > 0 ){
    sdf = as(sample_data(physeq), "data.frame")[, wh.svars, drop = FALSE]
    sdt = data.table(sdf, keep.rownames = TRUE)
    setnames(sdt, "rn", "SampleID")
    setkey(sdt, "SampleID")
    setkey(mdt, "SampleID")
    mdt <- sdt[mdt]
  }
  setkey(mdt, "TaxaID")
  return(mdt)
}

# Prevalences of each feature (stored as data frame)
prevdf <- apply(X = otu_table(ps),
               MARGIN = ifelse(taxa_are_rows(ps), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})

# Add taxonomy and total read counts to data frame
prevdf <- data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(ps),
                    tax_table(ps))

# Fast melt the phyloseq object
# Note: some downstream steps likely to be non-essential as they overlap with
# some of the clean-up steps in other tools. However, included here for completeness.

dt <- fast_melt(ps)

# Omit NAs
dt <- dt[!is.na(count)]

# Calculate relative abundances
dt[, RelativeAbundance := count / sum(count), by = SampleID]

# Make table with prevalence + total counts by OTU ID
prevdt <- dt[, list(Prevalence = sum(count > 0), TotalCounts = sum(count)), by = TaxaID]

# Count low-prevalence features
# (Prevalence values of 0-2) 
prev0 <- prevdt[(Prevalence <= 0), .N]
prev1 <- prevdt[(Prevalence <= 1), .N]
prev2 <- prevdt[(Prevalence <= 2), .N]

# Also count features with low total counts
# (Zero sequences, singletons, doubletons)

zeros <- prevdt[(TotalCounts <= 0), .N]
singletons <- prevdt[(TotalCounts <= 1), .N]
doubletons <- prevdt[(TotalCounts <= 2), .N]

# Subset prevdf for creating Phylum-level prevalence plot
prevdf_phyla <- subset(prevdf, Phylum %in% get_taxa_unique(ps, "Phylum"))

# Open a report PDF
pdf("ps_prevalence.pdf", width = 9, height = 7)

# Prevalence plot with 5 % filtering threshold included as guess
ggplot(prevdf_phyla, 
	aes(TotalAbundance, 
		Prevalence / nsamples(ps), 
		color = Phylum)) +
	geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) + 
	geom_point(size = 3, alpha = 0.7) +
	scale_x_log10() +
	facet_wrap(~ Phylum) + 
	theme(legend.position = "none") +
	theme(axis.text = element_text(size = 16)) +
	theme(axis.title = element_text(size = 20)) +
	theme(strip.text.x = element_text(size = 13)) +
	theme(axis.text.x = element_text(size = 16)) +
	theme(axis.text.x = element_text(size = 16, 
		angle = 45, vjust = 1, hjust = 1)) +
	xlab("\nTotal abundance") +
	ylab("Prevalence [frac. samples]\n")

# Close the report PDF
dev.off()

# Print out summary table
sink("ps_low.txt")
	cat("\n\n\n")
	cat("### Phyloseq summary: ###\n")
	cat("### Features with low prevalence and total counts ###\n")
	cat("\n\n\n")
	cat("### No.s of low-prevalence features ###\n")
	cat("\n\n\n")
	cat("### Prevalence = 0: ###\n")
	print(prev0)
	cat("\n\n\n")
	cat("### Prevalence = 1 (or less): ###\n")
	print(prev1)
	cat("\n\n\n")
	cat("### Prevalence = 2 (or less): ###\n")
	print(prev2)
	cat("\n\n\n")
	cat("### No.s of features with low total counts ###\n")
	cat("### (Including singletons and doubletons) ###\n")
	cat("\n\n\n")
	cat("### Zero sequences: ###\n")
	print(zeros)
	cat("\n\n\n")
	cat("### Total count = 1 (or less): ###\n")
	print(singletons)
	cat("\n\n\n")
	cat("### Total count = 2 (or less): ###\n")
	print(doubletons)
	cat("\n\n\n")
sink()
