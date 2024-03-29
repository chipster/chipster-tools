# TOOL ngs-find-nearest-genes.R: "Find the nearest genes for regions" (This tool takes a set of genomic regions, such as ChIP-seq peaks, and fetches the nearest gene for each.)
# INPUT regions-list.tsv: "Table with genomic regions" TYPE GENERIC
# OUTPUT nearest-genes.tsv: "Table listing the nearest gene feature for each input region."
# PARAMETER species: "Genome" TYPE [Human_hg18: "Human (hg18\)", Human_hg19: "Human (hg19\)", Mouse_mm9: "Mouse (mm9\)", Mouse_mm10: "Mouse (mm10\)", Rat_rn4: "Rat (rn4\)", Rat_rn5: "Rat (rn5\)", Zebrafish_Zv8: "Zebrafish (Zv8\)", Zebrafish_Zv9: "Zebrafish (Zv9\)"] DEFAULT Human_hg19 (The genome to use for fetching annotations.)
# PARAMETER chr_column: "Chromosome column" TYPE COLUMN_SEL DEFAULT chr (Column containing chromosome infomration of peaks)
# PARAMETER start_column: "Start coordinate column" TYPE COLUMN_SEL DEFAULT start (Column containing start coordinates of peaks)
# PARAMETER end_column: "End coordinate column" TYPE COLUMN_SEL DEFAULT end (Column containing end coordinates of peaks)

# 28.10.2010, MG Created
# 09.02.2012, EK
# 08.05.2014, MK New genomes added to the list. Added BED support.
# 11.09.2014, EK Clarified the parameters
# 01.08.2019, ML Fixed the order of column headers of the output file

# Possible additional parameters and output
# OUTPUT unique-genes-list.tsv: "Table listing the unique genes that can be mapped with gene symbols and entrez gene ids."
# PARAMETER feature: "Genomic feature" TYPE [TSS, exon, miRNA] DEFAULT TSS (The genomic feature to look for in each region. Currently, transcription start site, exon or miRNA are supported feature types)
# PARAMETER e.value.cutoff: "E-value cutoff" TYPE DECIMAL FROM 0 TO 100 DEFAULT 0.01 (This parameter controls the alignment stringency, where a lower value means better alignment.)

# Load the libraries
library(ChIPpeakAnno)
library(biomaRt)

if (species == "EMPTY") {
    stop("CHIPSTER-NOTE: Please choose genome that is to be used in the analysis.")
}

# Load the annotation data
if (species == "Human_hg18") {
    data(TSS.human.NCBI36)
    annotations <- TSS.human.NCBI36
    ensembl_dataset <- "hsapiens_gene_ensembl"
}
if (species == "Human_hg19") {
    data(TSS.human.GRCh37)
    annotations <- TSS.human.GRCh37
    ensembl_dataset <- "hsapiens_gene_ensembl"
}
if (species == "Mouse_mm9") {
    data(TSS.mouse.NCBIM37)
    annotations <- TSS.mouse.NCBIM37
    ensembl_dataset <- "mmusculus_gene_ensembl"
}
if (species == "Mouse_mm10") {
    data(TSS.mouse.GRCm38)
    annotations <- TSS.mouse.GRCm38
    ensembl_dataset <- "mmusculus_gene_ensembl"
}
if (species == "Rat_rn4") {
    data(TSS.rat.RGSC3.4)
    annotations <- TSS.rat.RGSC3.4
    ensembl_dataset <- "rnorvegicus_gene_ensembl"
}
if (species == "Rat_rn5") {
    data(TSS.rat.Rnor_5.0)
    annotations <- TSS.rat.Rnor_5.0
    ensembl_dataset <- "rnorvegicus_gene_ensembl"
}
if (species == "Zebrafish_Zv8") {
    data(TSS.zebrafish.Zv8)
    annotations <- TSS.zebrafish.Zv8
    ensembl_dataset <- "drerio_gene_ensembl"
}
if (species == "Zebrafish_Zv9") {
    data(TSS.zebrafish.Zv9)
    annotations <- TSS.zebrafish.Zv9
    ensembl_dataset <- "drerio_gene_ensembl"
}

# Read in data from BED or tsv file and convert it to BED format
if (length(grep("^column\\d+$", chr_column)) == 1 && length(grep("^column\\d+$", start_column)) == 1 && length(grep("^column\\d+$", end_column)) == 1) {
    results_file <- read.table(file = "regions-list.tsv", sep = "\t", header = F)
    chr_column <- as.numeric(gsub("^column", "", chr_column)) + 1
    start_column <- as.numeric(gsub("^column", "", start_column)) + 1
    end_column <- as.numeric(gsub("^column", "", end_column)) + 1
} else {
    results_file <- read.table(file = "regions-list.tsv", sep = "\t", header = T)
    chr_column <- grep(paste("^", chr_column, "$", sep = ""), colnames(results_file))
    start_column <- grep(paste("^", start_column, "$", sep = ""), colnames(results_file))
    end_column <- grep(paste("^", end_column, "$", sep = ""), colnames(results_file))
}

results_bed <- results_file[, c(chr_column, start_column, end_column)]
if (length(grep("chr", levels(results_bed[, 1]), invert = T)) > 0) {
    levels(results_bed[, 1])[grep("chr", levels(results_bed[, 1]), invert = T)] <- paste("chr", levels(results_bed[, 1])[grep("chr", levels(results_bed[, 1]), invert = T)], sep = "")
}

# Convert to Ranged data
results_ranged <- BED2RangedData(results_bed)

# Annotate the data
results_annotated <- annotatePeakInBatch(results_ranged, AnnotationData = annotations)

# Convert into a table and extract the useful columns
results_table <- as.data.frame(results_annotated)
# rownames(results_table) <- sub (" ", "_", results_table$name)
results_table <- results_table[, c(-5, -13, -14)]
table_header <- c("chromosome", "region_start", "region_end", "region_width", "region_id", "ensembl_id", "gene_start", "gene_end", "strand", "location", "distance")
names(results_table) <- table_header

# Write output
# write.table(results_table, file="nearest-genes.tsv", sep="\t", col.names=T, row.names=T, quote=F)
write.table(results_table, file = "nearest-genes.tsv", sep = "\t", col.names = T, row.names = F, quote = F)
# EOF
