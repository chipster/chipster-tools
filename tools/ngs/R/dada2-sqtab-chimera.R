# TOOL dada2-sqtab-chimera.R: "Make an ASV table and remove chimeras" (This tool makes an ASV-table and removes chimeras. As the input you can give either if you have single end reads the output .Rda object from the tool Sample inference or if you have paired end data the output .Rda object from the tool Combine paired reads to contigs with Dada2. )
# INPUT object.Rda: "Either mergers object named contigs.Rda or dada-class object named dada-forward.Rda"  TYPE GENERIC (Dada-class object named dada-forward.Rda if single end reads or a mergers object named contigs.Rda if paired end reads.)
# OUTPUT seqtab_nochim.Rda
# OUTPUT reads_summary.tsv
# OUTPUT summary.txt
# OUTPUT sequence_table_nochim.tsv
# PARAMETER method1: "Method to identify chimeras" TYPE [consensus, pooled] DEFAULT consensus (Identification by consensus across samples or identification from pooled sequences)
# RUNTIME R-4.4.3-asv
# TOOLS_BIN ""

source(file.path(chipster.common.lib.path, "tool-utils.R"))
source(file.path(chipster.common.lib.path, "zip-utils.R"))

# load library dada2
library(dada2)

load("object.Rda", verbose = TRUE)
# 2 possible names either mergers or dadaFs

if (exists("dadaFs")) { # single reads, rename
    object <- dadaFs
    name <- "After sample inference"
} else if (exists("mergers")) {
    object <- mergers
    name <- "After make contigs"
} else {
    stop(paste("CHIPSTER-NOTE: ", "The given .Rda object is not correct. If you are having single reads, use the output file from the tool Sample inference as the input and otherwise the output from the tool Combine paired reads to contigs with Dada2 "))
}
# Construct sequence table from given .Rda object -> ASV-table Construct a sample-by-sequence observation matrix.
seqtab <- makeSequenceTable(object)
# print(names(object))
# if mergers and one sample: "sequence"  "abundance" "forward"   "reverse"   "nmatch"    "nmismatch" "nindel"    "prefer"    "accept"
# if dadafs and one sample: "denoised"   "clustering" "sequence"   "quality" "birth_subs" "trans" "map" "err_in"  "err_out" "opts" "pval"
if (names(object)[1] != "sequence" && names(object)[1] != "denoised") {
    # rename the rows by splitting the name with _
    sample.names <- sapply(strsplit(rownames(seqtab), "_"), `[`, 1)
    rownames(seqtab) <- sample.names
}

# Distribution of sequence lengths, make it look nicer/ add rowname Counts:
data_table <- table(nchar(getSequences(seqtab)))
data <- t(data.frame(as.numeric(data_table))) # transpose of the vector
colnames(data) <- as.numeric(names(data_table)) # sequence lengths
rownames(data) <- "Counts:"

# Write a log/summary file
sink(file = "summary.txt")
cat("\nAfter  dada() algorithm sequence table consist of:\n")
cat(length(rownames(seqtab)), " samples and ", length(colnames(seqtab)), " amplicon sequence variants\n\n")
cat("Distribution of the amplicon sequence variant's lengths: Column names are the sequence lengths\n\n")
print(data)
cat("\n###Removing Chimeras:###\n")

# run isbimeradenovo / remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method = method1, multithread = as.integer(chipster.threads.max), verbose = TRUE)
num <- length(colnames(seqtab)) - length(colnames(seqtab.nochim))
cat("Identified ", num, " bimeras out of ", length(colnames(seqtab)), " input sequences\n")
cat("Total amount of ASVs is: ")
cat(length(colnames(seqtab.nochim)))
cat("\n\n")
# if (!is.na(mock)){
#  unqs.mock <- seqtab.nochim[mock,]
# unqs.mock <- sort(unqs.mock[unqs.mock>0], decreasing=TRUE) # Drop ASVs absent in the Mock
#  cat("DADA2 inferred", length(unqs.mock), "sample sequences present in the Mock community.\n")}
sink()


# track reads that were removed and make a tsv table
#
# rename sequences to asv1, asv2... makes visualization easier
getN <- function(x) sum(getUniques(x))
seqtab.nochim2 <- seqtab.nochim

if (names(object)[1] != "sequence" && names(object)[1] != "denoised") { # If processing a single sample, remove the sapply calls: e.g. replace sapply(object, getN) with getN(dadaFs)
    track <- cbind(sapply(object, getN), rowSums(seqtab.nochim))
    rownames(track) <- sample.names
    colnames(track) <- c(name, "Removed chimeras")
    write.table(track, file = "reads_summary.tsv", sep = "\t", row.names = TRUE, col.names = T, quote = F)

    colnames(seqtab.nochim2) <- paste0("ASV", seq(length(colnames(seqtab.nochim2))))
    write.table(seqtab.nochim2, file = "sequence_table_nochim.tsv", sep = "\t", row.names = TRUE, col.names = T, quote = F)
} else {
    track <- cbind(getN(object), rowSums(seqtab.nochim))
    # sample.names="Sample1" # rename to sample 1 if only one sample
    colnames(track) <- c(name, "Removed chimeras")
    write.table(track, file = "reads_summary.tsv", sep = "\t", row.names = FALSE, col.names = T, quote = F)

    colnames(seqtab.nochim2) <- paste0("ASV", seq(length(colnames(seqtab.nochim2))))
    write.table(seqtab.nochim2, file = "sequence_table_nochim.tsv", sep = "\t", row.names = FALSE, col.names = T, quote = F)
}




# print out asv sequence table rename asv:s
# rename sequences to asv1, asv2... makes it easier just for visualisation
# seqtab.nochim2 <- seqtab.nochim


# save the object as .Rda
save(seqtab.nochim, file = "seqtab_nochim.Rda")
