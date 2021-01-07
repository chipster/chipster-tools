# TOOL metabarcoding-transform-otu.R: "Transform OTU counts" (This tool is used to transform raw OTU counts within a phyloseq object. Four options are available\: 1\) CLR transformation \(with pseudocount\), 2\) relative abundance \(%\) conversion\; 3\) Hellinger transformation\; 4\) DESeq2 format conversion and variance-stabilizing transformation. The CLR transformation is based on the R package microbiome and applies a pseudocount of min\(relative abundance\)\/2 to zero relative abundance entries in the OTU table. Converting to DESeq2 format requires a user-specified phenodata variable as part of the DESeq2 experimental design formula. Note that data subjected to VST are saved as a phyloseq format \(data saved in the DESeq2 format are not VST-transformed\). Requires a phyloseq object in Rda format as the input. The resulting data are saved as an Rda file.)
# INPUT ps.Rda: "Phyloseq object in Rda format" TYPE GENERIC
# INPUT META phenodata.tsv: "Phenodata" TYPE GENERIC
# OUTPUT OPTIONAL ps_clr.Rda
# OUTPUT OPTIONAL ps_relabund.Rda
# OUTPUT OPTIONAL ps_hellinger.Rda
# OUTPUT OPTIONAL ps_vst.Rda
# OUTPUT OPTIONAL deseq2.Rda
# PARAMETER treatment: "Data treatment" TYPE [clr: "Centered log-ratio transfomation with pseudocount", relabund: "Relative abundances \(%\)", hellinger: "Hellinger transformation", deseq2: "DESeq2 format conversion and variance-stabilizing transformation"] DEFAULT clr (Choice between data transformation types)
# PARAMETER OPTIONAL group_column1: "Phenodata variable used for DESeq2 conversion" TYPE METACOLUMN_SEL DEFAULT empty (Select a phenodata variable used to specify the experimental design when converting the data to DESeq2 format.)
# RUNTIME R-3.6.1-phyloseq

# JH 2020

# Load packages
library(phyloseq)
library(microbiome)
library(DESeq2)

# Load phyloseq object
load("ps.Rda")

# CLR transformation

# Note:
# CLR transformation in microbiome package applies a pseudocount of min(relative abundance)/2 
# to exact zero relative abundance entries in OTU table before taking logs

if (treatment == "clr"){
	ps <- microbiome::transform(ps, 'clr')
}

# Relative abundance conversion
if (treatment == "relabund"){
	ps <- microbiome::transform(ps, 'compositional')
}

# Hellinger transformation
if (treatment == "hellinger"){
	ps <- microbiome::transform(ps, 'hellinger')
}

# DESeq2 conversion and variance-stabilizing transformation

# Note: size factors could be estimated e.g. by estimateSizeFactors(diagdds)
# However, here we use a zero-tolerant variant of geometric mean to avoid
# error messages with certain data sets (when there is a high prevalence
# of sparsely sampled OTUs; see https://github.com/joey711/phyloseq/issues/387)

# VST-transformed data assigned to ps (used for analyses other than DESeq2)
# Data without VST transformation assigned to ps0 (used for DESeq2 analysis)

if (treatment == "deseq2"){
	deseq2_formula <- as.formula(paste0("formula(~", group_column1, ")"))
	diagdds <- phyloseq_to_deseq2(ps, deseq2_formula)
	gm_mean <- function(x, na.rm = TRUE){
		exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))
		}
	geoMeans <- apply(counts(diagdds), 1, gm_mean)
	diagdds <- estimateSizeFactors(diagdds, geoMeans = geoMeans)
	diagdisp <- estimateDispersions(diagdds)
	diagvst <- getVarianceStabilizedData(diagdisp)
	ps0 <- ps
	otu_table(ps) <- otu_table(diagvst, taxa_are_rows = TRUE)
}

# Export CLR-transformed data
if (treatment == "clr"){
	save(ps, file = "ps_clr.Rda")
}

# Export relative abundance data
if (treatment == "relabund"){
	save(ps, file = "ps_relabund.Rda")
}

# Export Hellinger-transformed data
if (treatment == "hellinger"){
	save(ps, file = "ps_hellinger.Rda")
}

# Export VST-transformed phyloseq object;
# also export data in DESeq2 format (without VST)
if (treatment == "deseq2"){
	save(ps, file = "ps_vst.Rda")
	save(list = c("ps0", 
		"group_column1", "diagdds"), 
		file = "deseq2.Rda")
}
