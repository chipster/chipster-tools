# TOOL dada-make-phyloseq-object.R: "Make a Phyloseq object" (This tool makes a phyloseq object and saves it as a .Rda file. Requires an otu/asv-table which can be produced with the tool Make contigs and remove chimeras and a taxonomy table as the input.)
# INPUT chimera.Rda: "Saved class object produced with tool Make contigs and remove chimeras, named seqtab_nochim.Rda" TYPE GENERIC
# INPUT assignment.Rda: "Saved class object produced by assign taxonomy, named taxonomy-assignment-matrix.Rda" TYPE GENERIC
# OUTPUT META phenodata.tsv: "Phenodata" 
# OUTPUT ps_nophe.Rda
# OUTPUT ps_summary.txt
# PARAMETER OPTIONAL mock: "Name of the mock community if you have co-sequenced a mock community" TYPE STRING (If you have co-sequenced a mock community, you can remove it from the phyloseq object by giving the name of the community as a parameter. Most likely to be named as Mock)
# RUNTIME R-3.6.1-phyloseq

# RUNTIME R-4.1.1  no need for dada2 anymore


source(file.path(chipster.common.path,"tool-utils.R"))
source(file.path(chipster.common.path,"zip-utils.R"))

# Load phyloseq
library(phyloseq)
#packageVersion("phyloseq")
library(Biostrings)

# load input files
load("chimera.Rda", verbose=TRUE)
 # seqtab.nochim
load("assignment.Rda", verbose=TRUE)
 # taxa

# make a phyloseq object
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
      tax_table(taxa)) #sample_data(phenodata),

# Remove mock sample if it's selected
#if (!is.na(mock)){
#    ps <- prune_samples(sample_names(ps) != mock, ps)
#}

dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))


# make a summary text file
sink("ps_summary.txt")
	cat("\n\n\n")
	cat("### Phyloseq object ###\n")
	cat("\n\n\n")
	print(ps)
	cat("\n\n")
    #if (!is.na(mock)){
     #   cat("### Mock sequence:\n")
     #   unqs.mock <- seqtab.nochim[mock,]
     #   unqs.mock <- sort(unqs.mock[unqs.mock>0], decreasing=TRUE) # Drop ASVs absent in the Mock
     #   cat("DADA2 inferred", length(unqs.mock), "sample sequences present in the Mock community.\n")}
sink()

# Rename tax table columns to make the Phyloseq object suitable for OTU analyse tools
# Note 1: Rank1 renamed as Domain_Kingdom to reflect reference database-specific
# differences (some specify Rank1 as domain, others as kingdom)
# Note 2: Renaming depends on no. of taxonomic levels used (can be either 6 or 7 levels
taxlength <- length(colnames(tax_table(ps)))

if (taxlength == "6"){
colnames(tax_table(ps)) <- c("Domain_Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
}
if (taxlength == "7"){
colnames(tax_table(ps)) <- c("Domain_Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
}
# take out the sample names
samples.out <- sample_names(ps) 
# write phenodata.tsv
write.table(data.frame(sample = samples.out,
    chiptype = "NGS",
    group = rep("",length(samples.out))),
    "phenodata.tsv",
    col.names = T,
    row.names = F,
    sep = "\t",
    quote = F)



save(ps, file="ps_nophe.Rda")
#EOF