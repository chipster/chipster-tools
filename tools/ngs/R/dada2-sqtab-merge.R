# TOOL dada2-sqtab-merge.R: "Merge two ASV tables" (This tool merges two ASV tables, for example from different sequencing runs. To merge more than two ASV tables, run the tool several times.)
# INPUT seqtab_nochim1.Rda: "ASV table 1 saved as .Rda file." TYPE GENERIC (File is produced with tool "Make contigs and remove chimeras" and named seqtab_nochim.Rda, or it is a merged ASV table named seqtab_merged.Rda)
# INPUT seqtab_nochim2.Rda: "ASV table 2 saved as .Rda file." TYPE GENERIC (File is produced with tool "Make contigs and remove chimeras" and named seqtab_nochim.Rda, or it is a merged ASV table named seqtab_merged.Rda)
# OUTPUT seqtab_merged.Rda
# OUTPUT summary.txt
# OUTPUT sequence_table_merged.tsv
# RUNTIME R-4.4.3-asv
# TOOLS_BIN ""


# HJ 2025

source(file.path(chipster.common.lib.path, "tool-utils.R"))
source(file.path(chipster.common.lib.path, "zip-utils.R"))

# load library dada2
library(dada2)

# load and rename ASV tables
load("seqtab_nochim1.Rda", verbose = TRUE)
seqtab_nochim1 <- seqtab.nochim
load("seqtab_nochim2.Rda", verbose = TRUE)
seqtab_nochim2 <- seqtab.nochim
rm(seqtab.nochim)

# merge ASV tables
seqtab.nochim <- mergeSequenceTables(seqtab_nochim1, seqtab_nochim2)

# write a log/summary file
sink(file = "summary.txt")
cat("\nAfter merging, the sequence table consist of:\n")
cat(length(rownames(seqtab.nochim)), " samples and ", length(colnames(seqtab.nochim)), " amplicon sequence variants\n\n")
sink()

# save the object as .Rda
save(seqtab.nochim, file = "seqtab_merged.Rda")

# write the merged ASV table
colnames(seqtab.nochim) <- paste0("ASV", seq(length(colnames(seqtab.nochim))))
write.table(seqtab.nochim, file = "sequence_table_merged.tsv", sep = "\t", row.names = TRUE, col.names = T, quote = F)

