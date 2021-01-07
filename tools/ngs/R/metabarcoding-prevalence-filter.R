# TOOL metabarcoding-prevalence-filter.R: "Proportional prevalence filtering" (Filters out any OTUs below a user-specified prevalence threshold \(%\) using phyloseq. The resulting object is saved as an Rda file. Requires a phyloseq object in Rda format as the input.)
# INPUT ps.Rda: "Phyloseq object in Rda format" TYPE GENERIC
# OUTPUT ps_prevfilter.Rda: ps_prevfilter.Rda
# PARAMETER OPTIONAL threshold: "Filtering threshold (%)" TYPE INTEGER FROM 1 TO 25 DEFAULT 5 (Prevalence threshold \(%\) for filtering, accepts integer numbers from 1 to 25)
# RUNTIME R-3.6.1-phyloseq

# JH 2020

# Load phyloseq
library(phyloseq)

# Load phyloseq object
load("ps.Rda")

# Prevalences of each feature (stored as data frame)
prevdf <- apply(X = otu_table(ps),
               MARGIN = ifelse(taxa_are_rows(ps), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})

# Add taxonomy and total read counts to data frame
prevdf <- data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(ps),
                    tax_table(ps))

# Filter the data set
threshold <- threshold / 100
prevalenceThreshold <- threshold * nsamples(ps)
keepTaxa <- rownames(prevdf)[(prevdf$Prevalence >= prevalenceThreshold)]
ps <- prune_taxa(keepTaxa, ps)

# Export phyloseq object as Rda file 
save(ps, file = "ps_prevfilter.Rda")
