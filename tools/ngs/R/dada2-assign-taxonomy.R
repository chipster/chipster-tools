# TOOL dada2-assign-taxonomy.R: "Assign taxonomy" (Assign taxonomy to the sequence-variants.)
# INPUT seqtab_nochim.Rda: "Saved class object produced with tool Make contigs and remove chimeras, named seqtab_nochim.Rda" TYPE GENERIC
# INPUT OPTIONAL reference.fasta: "Own reference file" TYPE GENERIC
# OUTPUT taxonomy-assignment-matrix.Rda
# OUTPUT taxa_seqtab_combined.tsv
# OUTPUT taxonomy_assignment.tsv
# RUNTIME R-4.1.1

source(file.path(chipster.common.path,"tool-utils.R"))
source(file.path(chipster.common.path,"zip-utils.R"))

# load library dada2
library(dada2)
#packageVersion("dada2")

load("seqtab_nochim.Rda", verbose=TRUE)
#name seqtab.nochim
if (file.exists("reference.fasta")){
    path <- "reference.fasta"
}


# paths to silva-reference files
path <- c(file.path(chipster.tools.path,"dada2-silva-reference","silva_nr99_v138.1_train_set.fa"))
path2 <-  c(file.path(chipster.tools.path,"dada2-silva-reference","silva_species_assignment_v138.1.fa"))

# run command assignTaxonomy verbose not important
taxa <- assignTaxonomy(seqtab.nochim, path, multithread=TRUE)

# run addSpecies()
taxa <- addSpecies(taxa, path2)

# save the taxonomy table as taxonomy-assignment-matrix.Rda
save(taxa, file="taxonomy-assignment-matrix.Rda")

# rename seqtab.nochim colnames to asv1, asv2 instead of long sequences... 
names <- c()
x=0
while (x<length(colnames(seqtab.nochim))){
    new = paste("asv",x,sep="")
    names <- c(names, new)
    x = x+1
}
colnames(seqtab.nochim)<-names

# rename taxa rownames to asv1, asv2... 
names <- c()
x=0
while (x<length(rownames(taxa))){
    new = paste("asv",x,sep="")
    names <- c(names, new)
    x = x+1
}
rownames(taxa)<-names

# write taxa-table with rownames asv0...
write.table(taxa, file="taxonomy_assignment.tsv", sep="\t", row.names=TRUE, col.names=T, quote=F)

# combine taxa and seqtab.nochim matrixes to one table
df.combined <- cbind(taxa, t(seqtab.nochim))
write.table(df.combined, file="taxa_seqtab_combined.tsv", sep="\t", row.names=TRUE, col.names=T, quote=F)

#EOF