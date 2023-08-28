# TOOL dada-make-phyloseq-object.R: "Make a phyloseq object" (This tool makes a phyloseq object and saves it as a .Rda file. As the input, it requires an asv-table which can be produced with the tool /Make a sequence table and remove chimeras and a taxonomy table produced with the tool /Assign taxonomy)
# INPUT chimera.Rda: "Sequence table saved as a .Rda object, named seqtab_nochim.Rda" TYPE GENERIC
# INPUT assignment.Rda: "Taxonomy table saved as a .Rda object, named taxonomy-assignment-matrix.Rda " TYPE GENERIC
# OUTPUT META phenodata.tsv: "Phenodata"
# OUTPUT ps_nophe.Rda
# OUTPUT ps_summary.txt
# RUNTIME R-4.2.0-phyloseq


# ES 18.8.2022
# RUNTIME R-4.1.1-asv  no need for dada2
# PARAMETER OPTIONAL mock: "Name of the mock community if you have co-sequenced a mock community" TYPE STRING (If you have co-sequenced a mock community, you can remove it from the phyloseq object by giving the name of the community as a parameter. Most likely to be named as Mock)

source(file.path(chipster.common.lib.path, "tool-utils.R"))
source(file.path(chipster.common.path, "zip-utils.R"))

# Load phyloseq
library(phyloseq)
# packageVersion("phyloseq")
library(Biostrings)

# load input files
load("chimera.Rda", verbose = TRUE)
# seqtab.nochim
load("assignment.Rda", verbose = TRUE)
# taxa

# make a phyloseq object
ps <- phyloseq(
  otu_table(seqtab.nochim, taxa_are_rows = FALSE),
  tax_table(taxa)
) # sample_data(phenodata),

# Remove mock sample if it's selected
# if (!is.na(mock)){
#    ps <- prune_samples(sample_names(ps) != mock, ps)
# }
# Rename the names of the asv and save sequences to refseq
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
# if (!is.na(mock)){
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

if (taxlength == "6") {
  colnames(tax_table(ps)) <- c("Domain_Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
}
if (taxlength == "7") {
  colnames(tax_table(ps)) <- c("Domain_Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
}

# Replacing NAs in taxonomy table with the best available classification:
# extracting taxonomy table from the phyloseq object
ps.tax <- tax_table(ps)

# vector for abbreviations of taxonomic levels
taxabb <- c("ki", "ph", "cl", "or", "fa", "ge")

# for loop for replacing NAs
# when an entry does not have an NA, it is saved as tmptax
# if the next entry is NA, the previous tmptax is used in paste

for (i in 1:nrow(ps.tax)) {
  for (j in 1:ncol(ps.tax)) {
    if (is.na(ps.tax[i, j]) == TRUE) {
      ps.tax[i, j] <- paste(tmptax, taxabb[j], sep = "_")
    } else {
      tmptax <- ps.tax[i, j]
    }
  }
}
# returning the modified taxonomy table to the phyloseq object
tax_table(ps) <- ps.tax

# take out the sample names
samples.out <- sample_names(ps)
# write phenodata.tsv
write.table(
  data.frame(
    sample = samples.out,
    chiptype = "NGS",
    group = rep("", length(samples.out))
  ),
  "phenodata.tsv",
  col.names = T,
  row.names = F,
  sep = "\t",
  quote = F
)



save(ps, file = "ps_nophe.Rda")
# EOF
