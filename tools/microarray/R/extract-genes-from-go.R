# TOOL extract-genes-from-go.R: "Extract genes from GO term" (Fetches the genes that belong to a given GO term. Note that you have to give the ID of the GO term in a format like GO:0005515.)
# INPUT normalized.tsv: normalized.tsv TYPE GENE_EXPRS
# INPUT META phenodata.tsv: phenodata.tsv TYPE GENERIC
# OUTPUT extracted-from-GO.tsv: extracted-from-GO.tsv
# PARAMETER match.term: "Match term" TYPE STRING DEFAULT empty (String to search for.)

# Extract the rows from a table that contain probes that map genes annotated
# to a given gene ontology term
#
# MG 6.10.2010

# Parameter settings (default) for testing purposes
# match.term<-"GO:0005515"

# Reads the chiptype from phenodata table
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

# Loads the correct annotation library
lib2 <- sub(".db", "", lib)
library(package = lib, character.only = T)

# Loads the data
file <- c("normalized.tsv")
dat <- read.table(file, header = T, sep = "\t", row.names = 1)

# Extract the mapping info
lib3 <- paste(lib2, "GO2ALLPROBES", sep = "")
env <- get(lib3)
go.2.probes <- as.list(env)

# Extract all probes belonging to query
probes.list <- go.2.probes[match.term]
probes.go <- unlist(probes.list)

# Extract probes in query list
probes.query <- rownames(dat)

# Extract common probes
match.indices <- match(probes.go, probes.query, nomatch = 0)
match.indices <- unique(match.indices[match.indices > 0])

# Extract data
dat2 <- dat[match.indices, ]

# Writing out the extracted data
write.table(dat2, file = "extracted-from-GO.tsv", sep = "\t", row.names = T, col.names = T, quote = F)
