# TOOL dada-merge-phenodata.phyloseq.R: "Merge phenodata to the phyloseq object" (Merge a phenodata.tsv with the sample information to the phyloseq object. )
# INPUT ps_phe.Rda: "Phyloseq object without phenodata/sample information" TYPE GENERIC
# INPUT META phenodata.tsv: "Class object produced by assign taxonomy" TYPE GENERIC
# OUTPUT ps_imported.txt
# OUTPUT ps.Rda
# PARAMETER samplevar: "Phenodata variable with sequencing sample IDs" TYPE METACOLUMN_SEL DEFAULT EMPTY (Phenodata variable with unique IDs for each community profile.)
# RUNTIME R-3.6.1-phyloseq

source(file.path(chipster.common.path,"tool-utils.R"))
source(file.path(chipster.common.path,"zip-utils.R"))

# Load phyloseq
library(phyloseq)
#packageVersion("phyloseq")

# load input files
load("ps_phe.Rda", verbose=TRUE)

# read the phenodata.tsv meta file
phenodata <- read.delim("phenodata.tsv")
rownames(phenodata) <- phenodata[, samplevar[1]]

# merge the phenodata table with the sample information to the phyloseq object with merge_phyloseq() command
ps <- merge_phyloseq(ps, sample_data(phenodata)) #sample_data(phenodata),

# Print out basic descriptors
ps_samplevars <- sample_variables(ps)
ps_samplenames <- sample_names(ps)

# make a summary text file
sink("ps_imported.txt")
	cat("\n\n\n")
	cat("### Imported phyloseq object ###\n")
	cat("\n\n\n")
	print(ps)
	cat("\n\n\n")
	cat("### Sample names ###\n")
    cat("\n\n\n") 
    print(ps_samplenames)
    cat("\n\n\n")
    cat("### Sample variables ###\n")
    cat("\n\n\n")
    print(ps_samplevars)
    cat("\n\n\n")
sink()

save(ps, file="ps.Rda")

#print(sample_data(ps))
#EOF