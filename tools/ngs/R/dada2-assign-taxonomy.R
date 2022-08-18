# TOOL dada2-assign-taxonomy.R: "Assign taxonomy" (Assign taxonomy to the sequence variants. This tool uses SILVA v.138.1 reference fastas for assignment if own reference files are not provided. The names of the ASV sequences are changed to ASV1, ASV2... just for visualisation. Check the manual for more information.)
# INPUT seqtab_nochim.Rda: "Seqtab object saved as .Rda file." TYPE GENERIC (File is produced with tool "Make contigs and remove chimeras" and named seqtab_nochim.Rda.)
# INPUT OPTIONAL taxa_reference.fasta: "Own reference training fasta for assignTaxonomy" TYPE GENERIC (Own reference file for kingdom-genus level assignments. Otherwise use the SILVA v.138.1 reference file)
# INPUT OPTIONAL species_reference.fasta: "Own reference training fasta for assignSpecies" TYPE GENERIC (Own reference file for exact, 100% identity, species level assignment.)
# OUTPUT taxonomy-assignment-matrix.Rda
# OUTPUT OPTIONAL taxa_seqtab_combined.tsv
# OUTPUT OPTIONAL taxonomy_assignment.tsv
# PARAMETER boot: "The minimum bootstrap confidence for assigning a taxonomic level" TYPE INTEGER FROM 0 DEFAULT 50 (The minimum bootstrap confidence for assigning a taxonomic level)
# PARAMETER species: "Exact species level assignment?" TYPE [yes, no] DEFAULT yes (Do you want to assign the sequences to the species level if there is an exact match 100% identity between ASVs and sequenced reference strains?)
# PARAMETER combine_tables: "Combine the taxonomy and the sequence table" TYPE [yes,no] DEFAULT yes (If set to yes, it combines the taxonomy and the sequence/ASV table into one .tsv file, otherwise the tsv file consist only of the taxonomy table.)
# RUNTIME R-4.1.1-asv
source(file.path(chipster.common.path,"tool-utils.R"))
source(file.path(chipster.common.path,"zip-utils.R"))

# load library dada2
library(dada2)
#packageVersion("dada2")

load("seqtab_nochim.Rda", verbose=TRUE)
#name seqtab.nochim

# paths to silva-reference files: path1 for assign_taxonomy and path2 for addSpecies

if (file.exists("taxa_reference.fasta")){
    path1 <- "taxa_reference.fasta"
} else{
    path1 <- c(file.path(chipster.tools.path,"dada2-silva-reference","silva_nr99_v138.1_train_set.fa"))
    #path1 <- c(file.path(chipster.tools.path,"dada2-silva-reference","silva_nr99_v138.1_wSpecies_train_set.fa"))
}
set.seed(100) # Initialize random number generator for reproducibility
# run command assignTaxonomy verbose not important
taxa <- assignTaxonomy(seqtab.nochim, path1, minBoot=boot, multithread=FALSE)


#check if parameter species yes, otherwise skip species level assignment. Use reference file if selected otherwise Silva v.138.1
if (species=="yes"){
    if ((file.exists("species_reference.fasta"))){
        path2 <- "species_reference.fasta"
    }else{
        if (file.exists("taxa_reference.fasta")){
            stop(paste('CHIPSTER-NOTE: ',"You didn't give a reference file for exact Species level assignment, but you wanted to use your own 
            reference file for assignTaxonomy and selected the addSpecies parameter yes"))
        }else{
            path2 <-  c(file.path(chipster.tools.path,"dada2-silva-reference","silva_species_assignment_v138.1.fa"))
        }
    
    }
    # run addSpecies()
    taxa <- addSpecies(taxa, path2)
}

# save the taxonomy table as taxonomy-assignment-matrix.Rda
save(taxa, file="taxonomy-assignment-matrix.Rda")

# rename seqtab.nochim colnames to asv1, asv2 instead of long sequences... 
colnames(seqtab.nochim)<-paste0("ASV", seq(length(colnames(seqtab.nochim))))
# rename taxa rownames to asv1, asv2... 


rownames(taxa)<-paste0("ASV", seq(length(rownames(taxa))))
if (combine_tables=="yes"){
    # combine taxa and seqtab.nochim matrixes to one table
    df.combined <- cbind(taxa, t(seqtab.nochim))
    write.table(df.combined, file="taxa_seqtab_combined.tsv", sep="\t", row.names=TRUE, col.names=T, quote=F)   
}else{
    # print out only the taxonomy table
    write.table(taxa, file="taxonomy_assignment.tsv", sep="\t", row.names=TRUE, col.names=T, quote=F)
}





#EOF