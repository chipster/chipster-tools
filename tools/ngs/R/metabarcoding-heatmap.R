# TOOL metabarcoding-heatmap.R: "Heatmap of community profiles" (Visualizes the abundance data of the most abundant taxa at a selected taxonomic level as a heatmap. The order of samples is based on ordination.)
# INPUT ps.Rda: "Phyloseq object in Rda format" TYPE GENERIC
# INPUT META phenodata.tsv: "Phenodata" TYPE GENERIC
# OUTPUT ps_heatmap.pdf
# PARAMETER xvar: "Phenodata variable with sequencing sample IDs" TYPE METACOLUMN_SEL DEFAULT EMPTY (Phenodata variable with unique IDs for each community profile.)
# PARAMETER type: "Level of biological organization to display on the heatmap" TYPE [otu: "OTU/ASV", genus: "Genus", family: "Family", order: "Order", class: "Class", phylum: "Phylum"] DEFAULT otu (Default: OTU/ASV)
# PARAMETER OPTIONAL thresh: "Number of most abundant taxa to include at the selected taxonomic level" TYPE INTEGER FROM 1 TO 300 DEFAULT 50 (1-300, default 50)
# PARAMETER OPTIONAL colour: "Select a colour scheme" TYPE [dark: "Light on dark", light: "Dark on light", contrast: "Contrast emphasized"] DEFAULT dark
# PARAMETER OPTIONAL ordering: "Order samples by variable" TYPE METACOLUMN_SEL DEFAULT EMPTY (Samples are ordered by the alphabetical order of the chosen phenodata variable. Overrides the default similarity-based order.)
# RUNTIME R-4.4.3-phyloseq
# TOOLS_BIN ""

# HJ 2025

# Load packages
library(phyloseq)
library(ggplot2)

# Load phyloseq object
load("ps.Rda")

if (colour == "dark") { 
    lowcol <- "#000033"
    highcol <- "#66CCFF"
    navalue <- "black"
} else if (colour == "light") {
    lowcol <- "#66CCFF"
    highcol <- "#000033"
    navalue <- "white"
} else {
    #contrast
    lowcol <- "#FFFFCC"
    highcol <- "#000033"
    navalue <- "white"
}

if (type == "otu") {
    ps <- prune_taxa(names(sort(taxa_sums(ps),TRUE)[1:thresh]), ps)
} else {

# Agglomerate data at the desired taxonomic level
# Functions also include some (likely superfluous) data cleaning
# Keep 'thresh' number of most abundant taxa
if (type == "phylum") {
    ps <- tax_glom(ps,
        taxrank = rank_names(ps)[2],
        NArm = TRUE, bad_empty = c(NA, "", " ", "\t"))
    taxa_names(ps) <- tax_table(ps)[,2]
    ps <- prune_taxa(names(sort(taxa_sums(ps),TRUE)[1:thresh]), ps)

} else if (type == "class") {
    ps <- tax_glom(ps,
        taxrank = rank_names(ps)[3],
        NArm = TRUE, bad_empty = c(NA, "", " ", "\t"))
    taxa_names(ps) <- tax_table(ps)[,3]
    ps <- prune_taxa(names(sort(taxa_sums(ps),TRUE)[1:thresh]), ps)

} else if (type == "order") {
    ps <- tax_glom(ps,
        taxrank = rank_names(ps)[4],
        NArm = TRUE, bad_empty = c(NA, "", " ", "\t"))
    taxa_names(ps) <- tax_table(ps)[,4]
    ps <- prune_taxa(names(sort(taxa_sums(ps),TRUE)[1:thresh]), ps)

} else if (type == "family") {
    ps <- tax_glom(ps,
        taxrank = rank_names(ps)[5],
        NArm = TRUE, bad_empty = c(NA, "", " ", "\t"))
    taxa_names(ps) <- tax_table(ps)[,5]
    ps <- prune_taxa(names(sort(taxa_sums(ps),TRUE)[1:thresh]), ps)

} else if (type == "genus") {
    ps <- tax_glom(ps,
        taxrank = rank_names(ps)[6],
        NArm = TRUE, bad_empty = c(NA, "", " ", "\t"))
    taxa_names(ps) <- tax_table(ps)[,6]
    ps <- prune_taxa(names(sort(taxa_sums(ps),TRUE)[1:thresh]), ps)

}
}

if (ordering == "EMPTY") { 
psheatmap <- plot_heatmap(ps, "NMDS", "bray", low = lowcol, high = highcol, na.value = navalue)
} else
psheatmap <- plot_heatmap(ps, "NMDS", "bray", low = lowcol, high = highcol, na.value = navalue,
    sample.order = ordering)

# Open a report PDF
pdf("ps_heatmap.pdf", width = 12, height = 7)

psheatmap

# Close the report PDF
dev.off()
