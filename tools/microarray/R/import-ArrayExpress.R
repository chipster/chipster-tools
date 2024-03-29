# TOOL import-ArrayExpress.R: "Import from ArrayExpress" (Import Affymetrix raw data directly from ArrayExpress, and normalise it using the RMA algorithm.)
# OUTPUT normalized.tsv: normalized.tsv
# OUTPUT META phenodata.tsv: phenodata.tsv
# PARAMETER accession: Accession TYPE STRING DEFAULT E-MEXP-1422 (Accession number of the experiment.)
# PARAMETER platform: Platform TYPE STRING DEFAULT EMPTY (In case the series contains multiple platforms, specify the accession of the platform to import. If there is just one, this platform is ignored.)

# JTT: 13.1.2010
# MG: 30.4.2010, add detection calls and remove offending characters
# MK: 01.11.2013, add possibility to choose platform

# Parameter settings (default) for testing purposes
# accession<-"E-MEXP-1422"

# Loads the libraries
library(ArrayExpress)

# Loads the data
dat <- ArrayExpress(accession)

# if dataset contains multiple affy-sets, class is list
if (class(dat) == "list") {
   list_index <- NULL
   for (i in 1:length(dat)) {
      if (annotation(dat[[i]]) == platform) {
         list_index <- i
      }
   }
   if (is.null(list_index)) {
      stop("CHIPSTER-NOTE: Platforms matching your query were not found from this dataset")
   }
   dat <- dat[[list_index]]
}

# Normalizes the raw data
if (class(dat) == "AffyBatch") {
   library(affy)
   dat2 <- rma(dat)
   dat2 <- exprs(dat2)
   calls <- exprs(mas5calls(dat))
   calls <- as.data.frame(calls)
   colnames(dat2) <- paste("chip.", colnames(dat2), sep = "")
   names(calls) <- paste("flag.", names(calls), sep = "")
   dat2 <- data.frame(dat2, calls)
   chiptype <- dat@annotation
}

# Writing out data
a <- try(library(paste(chiptype, ".db", sep = ""), character.only = T))
if (chiptype != "empty" & class(a) != "try-error") {
   # Including gene names to data
   lib2 <- sub(".db", "", chiptype)
   symbol <- gsub("\'", "", data.frame(unlist(as.list(get(paste(lib2, "SYMBOL", sep = "")))))[rownames(dat2), ])
   genename <- gsub("\'", "", data.frame(unlist(as.list(get(paste(lib2, "GENENAME", sep = "")))))[rownames(dat2), ])
   # Writes the results into a file
   write.table(data.frame(symbol, description = genename, dat2), file = "normalized.tsv", col.names = T, quote = F, sep = "\t", row.names = T)
}

if (chiptype == "empty" | class(a) == "try-error") {
   write.table(data.frame(dat2), file = "normalized.tsv", col.names = T, quote = F, sep = "\t", row.names = T)
}

# Writes out a phenodata table
sample <- rownames(pData(dat))
# remove any comment character from phenodata info
phenodata.info <- pData(dat)
for (count in 1:dim(phenodata.info)[2]) {
   phenodata.info[, count] <- gsub("#", "", phenodata.info[, count])
}

group <- c(rep("", nrow(pData(dat))))
chiptype <- paste(chiptype, ".db", sep = "")
write.table(data.frame(sample = sample, chiptype = chiptype, group = group, phenodata.info), file = "phenodata.tsv", sep = "\t", row.names = F, col.names = T, quote = F)
