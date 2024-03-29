# TOOL add-locations-to-data.R: "Add genomic location to data" (Annotates genes with chromosomal location information, and adds the results to the datafile. )
# INPUT normalized.tsv: normalized.tsv TYPE GENE_EXPRS
# INPUT META phenodata.tsv: phenodata.tsv TYPE GENERIC
# OUTPUT data-with-locations.tsv: data-with-locations.tsv
# PARAMETER annotate_with: "Annotate using" TYPE [probe_id: "probe ID", gene_symbol: "gene symbols"] DEFAULT gene_symbol (Should the probe identifiers be used to fetch the location of the corresponding gene targets or should gene symbols be used directly, if available.)
# PARAMETER species: Species TYPE [human: human, mouse: mouse, rat: rat] DEFAULT human (The species needs to be specified in order to map genes to the genomic coordinates.)

# MG 16.03.2012
# IS  9.6.2013 Fixed sort order and changed output column name from 'chr' to 'chromosome' to be compatible with copy number scripts.
# MK 10.06.2013, fixing biomaRt queries
# ML 19.10.2016, fixing

# Loads libraries into memory
library(biomaRt)

# Reads the chiptype from phenodata table
if (annotate_with == "probe_id") {
    phenodata <- read.table("phenodata.tsv", header = T, sep = "\t")
    if (phenodata$chiptype[1] != "cDNA" | phenodata$chiptype[1] != "Illumina") {
        # Saves the chiptype into object lib
        lib <- phenodata$chiptype[1]
        lib <- as.character(lib)
    }
    # Account for the fact that annotation packages are from version 2.3 of Bioconductor
    # named with an ".db" suffix. Add the suffix when missing to support data files
    # from Chipster 1.3 and earlier.
    if (length(grep(".db", lib)) == 0 & length(grep("pmcdf", lib)) == 0) {
        lib <- paste(lib, ".db", sep = "")
    }
    # Checks that the package actually exists
    is_package_available <- try(library(package = lib, character.only = T))
    if (length(grep("no package", is_package_available[1])) > 0) {
        stop("CHIPSTER-NOTE: It appears that your data was not acquired using one of the array types supported in Chipster. If the input data already has a 'symbol' column, please rerun the tool with the 'annotate using' parameter set to 'gene symbols'.")
    }
    # Loads the correct annotation library
    library(package = lib, character.only = T)
    library(annaffy)
}

# Loads the data
file <- c("normalized.tsv")
dat <- read.table(file, header = T, sep = "\t", row.names = 1)

# Do sanity check if gene symbol is to be used
if (annotate_with == "gene_symbol") {
    gene_symbols <- dat[, grep("symbol", names(dat))]
    if (length(unique(gene_symbols)) < 1) {
        stop("CHIPSTER-NOTE: You need to have at least one gene symbol in your data table to use that information for fetching chromosome location information. For supported array designs the user can opt to use the probe ID instead.")
    }
}

# Covert probe id to gene symbols using appropriate annotation package
probes_query <- rownames(dat)
if (annotate_with == "probe_id") {
    probes_query <- rownames(dat)
    lib2 <- sub(".db", "", lib)
    lib3 <- paste(lib2, "SYMBOL", sep = "")
    env <- get(lib3)
    gene_symbols <- unlist(mget(probes_query, env, ifnotfound = NA))
}

# Fetch the gene symbols and descriptions from ENSEMBL using biomaRt
if (species == "human") {
    dataset <- "hsapiens_gene_ensembl"
    filt <- "hgnc_symbol"
}
if (species == "mouse") {
    dataset <- "mmusculus_gene_ensembl"
    filt <- "mgi_symbol"
}
if (species == "rat") {
    dataset <- "rnorvegicus_gene_ensembl"
    filt <- "rgd_symbol"
}

ensembl <- useMart("ensembl", dataset = dataset)
annotated_genes <- getBM(mart = ensembl, attributes = c(filt, "chromosome_name", "start_position", "end_position"), filters = filt, values = gene_symbols)

# Remove chromosome entries containing the "_" character
annotated_genes <- annotated_genes[-grep(pattern = "_", annotated_genes$chromosome_name), ]

# Match the list of input gene ids with the annotations
chr_name <- character(length(probes_query))
chr_start <- character(length(probes_query))
chr_end <- character(length(probes_query))
for (gene_count in 1:length(gene_symbols)) {
    chr_name[gene_count] <- annotated_genes[annotated_genes$hgnc_symbol == gene_symbols[gene_count], 2][1]
    chr_start[gene_count] <- annotated_genes[annotated_genes$hgnc_symbol == gene_symbols[gene_count], 3][1]
    chr_end[gene_count] <- annotated_genes[annotated_genes$hgnc_symbol == gene_symbols[gene_count], 4][1]
}
result_table <- data.frame(chromosome = chr_name, start = chr_start, end = chr_end, dat, stringsAsFactors = FALSE)

# Order the output based on chromosome and start position
result_table$chromosome <- factor(result_table$chromosome, levels = c(1:22, "X", "Y"), ordered = TRUE)
result_table <- result_table[order(result_table$chromosome, result_table$start), ]


# Output the data table with added chromosome location info
write.table(result_table, file = "data-with-locations.tsv", sep = "\t", row.names = T, col.names = T, quote = F)

# EOF
