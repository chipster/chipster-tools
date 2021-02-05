# TOOL metabarcoding-rarecurve-alpha.R: "Sequence numbers, rarefaction curve and alpha diversity estimates" (Lists per-sample sequence numbers, plots rarefaction curves and tabulates alpha diversity estimates \(observed no. of OTUs, Chao1 and Shannon indices, and Pielou's evenness\). Requires a phyloseq object in Rda format as the input. Note that the diversity estimates are only reliable when using raw \(untrimmed\) OTU data, as many diversity metrics are dependent on singletons and doubletons in the data set under analysis.)
# INPUT ps.Rda: "Phyloseq object in Rda format" TYPE GENERIC
# INPUT META phenodata.tsv: "Phenodata" TYPE GENERIC
# OUTPUT ps_rarecurve.pdf
# OUTPUT ps_alphadiv.txt
# PARAMETER OPTIONAL group_column: "Phenodata variable tabulated with alpha diversity estimates" TYPE METACOLUMN_SEL DEFAULT empty (Phenodata variable added to alpha diversity table for improved readability.)
# PARAMETER OPTIONAL type: "Are the data raw (untrimmed) or rarefied to even depth?" TYPE [raw, rarefied] DEFAULT raw (Type of data to be analyzed \(raw vs rarefied; by default, assumes raw data.\))
# RUNTIME R-3.6.1-phyloseq

# JH 2020

# Load libraries
library(microbiome)
library(phyloseq)
library(vegan)

# Load data
load("ps.Rda")

# Calculate per-sample sequence no.s
seqno <- sample_sums(ps)

# Alpha diversity estimates and Pielou's evenness
if (group_column == "empty") {
	set.seed(1)
	richness <- estimate_richness(ps, measures = c("Observed", "Chao1", "Shannon"))
	set.seed(1)
	pielou <- evenness(ps, 'pielou')
	richness <- cbind(richness, pielou)
}
if (group_column != "empty") {
	set.seed(1)
	richness <- estimate_richness(ps, measures = c("Observed", "Chao1", "Shannon"))
	set.seed(1)
	pielou <- evenness(ps, 'pielou')
	richness <- cbind(richness, pielou)	
	colno <- pmatch(group_column, colnames(ps@sam_data)) # Get desired column no. from ps metadata
	richness <- cbind(richness, ps@sam_data[,colno]) # cbind the column to the diversity table
}

# Open a report PDF
pdf("ps_rarecurve.pdf")

# Plot rarefaction curve
set.seed(1)
rarecurve(t(otu_table(ps)), step = 100, 
          cex.lab = 1.5, cex.axis = 1.5, label = FALSE, ylab = "OTUs", xlab = "No. of sequences")

# Close the report PDF
dev.off()

# Print out sequence numbers and alpha diversity table
if (type == "raw"){
sink("ps_alphadiv.txt")
	cat("\n\n\n")
	cat("### Per-sample sequence no.s ###\n")
	cat("\n\n\n")
	print(seqno)
	cat("\n\n\n")
	cat("### Alpha diversity estimates (observed OTUs, Chao1, Shannon's index, Pielou's evenness) ###\n")
	cat("\n\n\n")
	print(richness)
	cat("\n\n\n")
sink()
}
if (type == "rarefied"){
sink("ps_alphadiv.txt")
	cat("\n\n\n")
	cat("### Per-sample sequence no.s ###\n")
	cat("\n\n\n")
	print(seqno)
	cat("\n\n\n")
	cat("### Alpha diversity estimates (observed OTUs, Chao1, Shannon's index, Pielou's evenness) ###\n")
	cat("\n\n\n")
	cat("N/A (unavailable due to use of rarefied data)\n")
	cat("\n\n\n")
sink()
}
