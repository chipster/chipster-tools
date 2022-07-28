# TOOL dada2-assign-taxonomy.R: "Assign taxonomy" (Assign taxonomy to the sequence variants. This tool uses SILVA v.138.1 reference fastas for assignment if own reference files are not provided. Check the manual for more information.)
# INPUT seqtab_nochim.Rda: "Seqtab object saved as .Rda file." TYPE GENERIC (File is produced with tool "Make contigs and remove chimeras" and named seqtab_nochim.Rda.)
# INPUT OPTIONAL taxa_reference.fasta: "Own reference training fasta for assignTaxonomy" TYPE GENERIC (Own reference file for kingdom-genus level assignments.)
# INPUT OPTIONAL species_reference.fasta: "Own reference training fasta for assignSpecies" TYPE GENERIC (Own reference file for species level assignments. )
# OUTPUT taxonomy-assignment-matrix.Rda
# OUTPUT taxa_seqtab_combined.tsv
# OUTPUT taxonomy_assignment.tsv
# PARAMETER species: "Species level assignment?" TYPE [yes, no] DEFAULT yes (Do you want to assign the sequences to the species level?)
# RUNTIME R-4.1.1

source(file.path(chipster.common.path,"tool-utils.R"))
source(file.path(chipster.common.path,"zip-utils.R"))

# load library dada2
library(dada2)
#packageVersion("dada2")

load("seqtab_nochim.Rda", verbose=TRUE)
#name seqtab.nochim
# paths to silva-reference files path1 for assign_taxonomy and path2 for addSpecies

if (file.exists("taxa_reference.fasta")){
    path1 <- "taxa_reference.fasta"
} else{
    path1 <- c(file.path(chipster.tools.path,"dada2-silva-reference","silva_nr99_v138.1_train_set.fa"))
}
# run command assignTaxonomy verbose not important
taxa <- assignTaxonomy(seqtab.nochim, path1, multithread=TRUE)


#check if parameter species yes, otherwise skip species level assignment. Use reference file if selected otherwise Silva v.138.1
if (species=="yes"){
    if ((file.exists("species_reference.fasta"))){
        path2 <- "species_reference.fasta"
    }else{
        path2 <- path2 <-  c(file.path(chipster.tools.path,"dada2-silva-reference","silva_species_assignment_v138.1.fa"))
    }
    # run addSpecies()
    taxa <- addSpecies(taxa, path2)
}

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