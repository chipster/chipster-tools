# TOOL phyloseq-info.R: "Extract information from the Phyloseq object" (You can use this tool to extract information from the Phyloseq object: OTU/ASV-table, taxonomy table, sample information, refseg stored DNA sequences)
# INPUT ps.Rda: "Phyloseq object" TYPE GENERIC
# OUTPUT OPTIONAL otu_table.tsv
# OUTPUT OPTIONAL taxonomy_table.tsv
# OUTPUT OPTIONAL dna_sequences.fasta
# OUTPUT OPTIONAL dna_sequences.tsv
# OUTPUT OPTIONAL sample_information.tsv
# PARAMETER otu: "Extract OTU table" TYPE [yes,no] DEFAULT no
# PARAMETER taxa: "Extract taxonomy table" TYPE [yes,no] DEFAULT no
# PARAMETER fasta: "Extract DNA sequences in FASTA" TYPE [yes,no] DEFAULT no
# PARAMETER ref: "Extract DNA sequences as tabular data" TYPE [yes,no] DEFAULT no
# PARAMETER sample: "Extract sample information" TYPE [yes,no] DEFAULT no
# RUNTIME R-4.2.0-phyloseq


# ES 15.8.2022 phyloseq from 1.30 -> 1.40
# HJ 25.2.2025 add export to fasta
# Load phyloseq

library(phyloseq)
# packageVersion("phyloseq")

if (otu == "no" && sample == "no" && taxa == "no" && ref == "no" && ref == "no" && fasta == "no") {
    stop(paste("CHIPSTER-NOTE: ", "You didn't select any parameter. You have to set at least one parameter to yes."))
}
# load input files
load("ps.Rda")
if (otu == "yes") {
    otu_table <- otu_table(ps)
    write.table(otu_table, file = "otu_table.tsv", sep = "\t", row.names = TRUE, col.names = T, quote = F)
}
if (taxa == "yes") {
    tax_table <- tax_table(ps)
    write.table(tax_table, file = "taxonomy_table.tsv", sep = "\t", row.names = TRUE, col.names = T, quote = F)
}
if (fasta == "yes") {
    if (is.null(refseq(ps, errorIfNULL = FALSE))) {
        stop(paste("CHIPSTER-NOTE: ", "The refseq object is empty"))
    } else {
        refseq(ps) |>
                Biostrings::writeXStringSet("dna_sequences.fasta", append=FALSE,
                        compress=FALSE, compression_level=NA, format="fasta")
    }
}
if (ref == "yes") {
    if (is.null(refseq(ps, errorIfNULL = FALSE))) {
        stop(paste("CHIPSTER-NOTE: ", "The refseq object is empty"))
    } else {
        refseq <- refseq(ps)
        print(refseq)
        write.table(refseq, file = "dna_sequences.tsv", sep = "\t", row.names = TRUE, col.names = "DNA Sequences", quote = F)
    }
}
if (sample == "yes") {
    if (is.null(sample_data(ps, errorIfNULL = FALSE))) {
        stop(paste("CHIPSTER-NOTE: ", "The sample_data object is empty"))
    } else {
        sample_data <- sample_data(ps)
        write.table(sample_data, file = "sample_information.tsv", sep = "\t", row.names = F, col.names = T, quote = F)
    }
}

# EOF
