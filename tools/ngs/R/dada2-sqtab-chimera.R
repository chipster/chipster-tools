# TOOL dada2-sqtab-chimera.R: "Make a sequence table and remove chimeras" (This tool makes a sequence table / ASV-table and removes chimeras. As the input you can give either if you have single end reads the output .Rda object from the tool Sample inference or if you have paired end data the output .Rda object from the tool Combine paired reads to contigs with Dada2. )
# INPUT object.Rda: "Either mergers object named contigs.Rda or dada-class object named dada-forward.Rda"  TYPE GENERIC (Dada-class object named dada-forward.Rda if single end reads and a mergers object named contigs.Rda if paired end reads.)
# OUTPUT seqtab_nochim.Rda
# OUTPUT reads_summary.tsv
# OUTPUT summary.txt
# OUTPUT sequence_table_nochim.tsv
# PARAMETER method1: "Method to identify chimeras" TYPE [consensus, pooled] DEFAULT consensus (Identification by consensus across samples or identification from pooled sequences)
# RUNTIME R-4.1.1


source(file.path(chipster.common.path,"tool-utils.R"))
source(file.path(chipster.common.path,"zip-utils.R"))

# load library dada2
library(dada2)

load("object.Rda",verbose=TRUE)
# 2 possible names either mergers or dadaFs

if(exists("dadaFs")){ #single reads, rename
    object <- dadaFs
    name <- "After sample inference"
}else if(exists("mergers")){
    object <- mergers
    name <- "After make contigs"
}else{
    stop(paste('CHIPSTER-NOTE: ',"The given .Rda object is not correct. If you are having single reads, use the output file from the tool Sample inference as the input and otherwise the output from the tool Combine paired reads to contigs with Dada2 "))
}
# Construct sequence table from given .Rda object -> ASV-table Construct a sample-by-sequence observation matrix.
seqtab <- makeSequenceTable(object)

#rename the rows by splitting the name with _
sample.names <- sapply(strsplit(rownames(seqtab), "_"), `[`, 1)
rownames(seqtab) <- sample.names

# Write a log/summary file
sink(file="summary.txt")
    cat("\nAfter  dada() algorithm sequence table consist of:\n")
    cat(length(rownames(seqtab)))
    cat(" samples and ")
    cat(length(colnames(seqtab)))
    cat(" amplicon sequence variants\n")
    cat("\nDistribution of sequence lengths:\n")
    table(nchar(getSequences(seqtab)))
    cat("\n###Removing Chimeras:###\n")

# run isbimeradenovo / remove chimeras
    seqtab.nochim <- removeBimeraDenovo(seqtab, method=method1, multithread=TRUE, verbose=TRUE)
    num <- length(colnames(seqtab))-length(colnames(seqtab.nochim))
    cat("Identified ",num," bimeras out of ",length(colnames(seqtab))," input sequences\n")
    cat("Total amount of ASVs is: ")
    cat(length(colnames(seqtab.nochim)))
    cat("\n\n")
    #if (!is.na(mock)){
    #  unqs.mock <- seqtab.nochim[mock,]
     # unqs.mock <- sort(unqs.mock[unqs.mock>0], decreasing=TRUE) # Drop ASVs absent in the Mock
    #  cat("DADA2 inferred", length(unqs.mock), "sample sequences present in the Mock community.\n")}
sink()

# track readsthat were removed asand make a tsv table
getN <- function(x) sum(getUniques(x))
track <- cbind(sapply(object, getN), rowSums(seqtab.nochim))
colnames(track) <- c( name, "Removed chimeras")
rownames(track) <- sample.names
write.table(track, file="reads_summary.tsv", sep="\t", row.names=TRUE, col.names=T, quote=F)


# print out asv sequence table rename asv:s
# rename sequences to asv1, asv2... makes it easier just for visualisation
seqtab.nochim2 <- seqtab.nochim
colnames(seqtab.nochim2)<-paste0("ASV", seq(length(colnames(seqtab.nochim2))))
write.table(seqtab.nochim2, file="sequence_table_nochim.tsv", sep="\t", row.names=TRUE, col.names=T, quote=F)

# save the object as .Rda 
save(seqtab.nochim, file = "seqtab_nochim.Rda")


print(seqtab.nochim[rownames(seqtab.nochim) != "Mock"],)