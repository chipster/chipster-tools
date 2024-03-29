# TOOL mothur-phyloseq-inputs.R: "Cluster sequences to OTUs and classify them" (This tool generates input files for the R package phyloseq, using the Mothur tools dist.seqs \(for 16S and 18S data\), cluster \(opticlust for 16S and 18S data, and agc for ITS data\), make.shared and classify.otu. As part of input file generation, sequences are clustered into operational taxonomic units \(OTU\) based on a user-defined sequence dissimilarity threshold. Note that sample groups in the phenodata file should be indicated using characters \(numeric labels are currently not supported by downstream analyses\).)
# INPUT file.fasta: "FASTA file" TYPE GENERIC
# INPUT picked.count_table: "Mothur count file" TYPE MOTHUR_COUNT
# INPUT sequences-taxonomy-assignment.txt: "Sequences taxonomy assignment file" TYPE GENERIC
# OUTPUT META phenodata.tsv
# OUTPUT OPTIONAL file.opti_mcc.0.05.cons.taxonomy
# OUTPUT OPTIONAL file.opti_mcc.0.04.cons.taxonomy
# OUTPUT OPTIONAL file.opti_mcc.0.03.cons.taxonomy
# OUTPUT OPTIONAL file.opti_mcc.0.02.cons.taxonomy
# OUTPUT OPTIONAL file.opti_mcc.0.01.cons.taxonomy
# OUTPUT OPTIONAL file.agc.0.05.cons.taxonomy
# OUTPUT OPTIONAL file.agc.0.04.cons.taxonomy
# OUTPUT OPTIONAL file.agc.0.03.cons.taxonomy
# OUTPUT OPTIONAL file.agc.0.02.cons.taxonomy
# OUTPUT OPTIONAL file.agc.0.01.cons.taxonomy
# OUTPUT OPTIONAL file.opti_mcc.shared
# OUTPUT OPTIONAL file.agc.shared
# PARAMETER datatype: "Type of data" TYPE [other: "16S or 18S", its: "ITS"] DEFAULT other (Indicate if you have ITS data as it is treated differently.)
# PARAMETER cutoff: "Clustering cutoff" TYPE [0.05, 0.04, 0.03, 0.02, 0.01] DEFAULT 0.03 (Dissimilarity threshold for OTU clustering, e.g. a cut-off value of 0.03 corresponds to 97% similarity)
# RUNTIME R-4.2.0-phyloseq
# SLOTS 4
# STORAGE 2000

# AMS JH EK 2020-2021
# ES 28.12.2022 not possible to update to new mothur version because reshape library not in the RUNTIME R-4.1.1, use old until reshsape updated
# EK 10.01.2023 updated to Mothur 1.48.0 using RUNTIME R-4.2.0-phyloseq which contains the reshape library

# reshape2 library
library(reshape2)

# check out if the file is compressed and if so unzip it
source(file.path(chipster.common.lib.path, "tool-utils.R"))
source(file.path(chipster.common.lib.path, "zip-utils.R"))
unzipIfGZipFile("file.fasta")

# binary
binary <- c(file.path(chipster.tools.path, "mothur", "mothur"))
version <- system(paste(binary, "--version"), intern = TRUE)
documentVersion("Mothur", version)

# mothur dist.seqs
# used for 16S and 18S data only

if (datatype == "other") {
    distseqs.options <- paste("dist.seqs(fasta=file.fasta") # dist.seqs produces file.dist
    distseqs.options <- paste(distseqs.options, ", processors=", chipster.threads.max, sep = "")
    distseqs.options <- paste(distseqs.options, ", cutoff=", cutoff, ")", sep = "")
    documentCommand(distseqs.options)
    write(distseqs.options, "distseqs.mth", append = FALSE)
    command <- paste(binary, "distseqs.mth", "> log_distseqs.txt")
    # system(command)
    runExternal(command, checkexit = TRUE)
}

# mothur cluster

# method depends on whether data are ITS or not
# (ITS contig lengths are not equal, whereas lengths are equal for other data)

if (datatype == "other") {
    cluster.options <- paste("cluster(column=file.dist, count=picked.count_table")
    cluster.options <- paste(cluster.options, ", cutoff=", cutoff, ")", sep = "", "seed" = 1)
    documentCommand(cluster.options)
    write(cluster.options, "cluster.mth", append = FALSE)
    command <- paste(binary, "cluster.mth", "> log_cluster.txt")
    # system(command)
    runExternal(command, checkexit = TRUE)
}

if (datatype == "its") {
    cluster.options <- paste("cluster(fasta=file.fasta, count=picked.count_table, method=agc")
    cluster.options <- paste(cluster.options, ", cutoff=", cutoff, ")", sep = "", "seed" = 1)
    documentCommand(cluster.options)
    write(cluster.options, "cluster.mth", append = FALSE)
    command <- paste(binary, "cluster.mth", "> log_cluster.txt")
    # system(command)
    runExternal(command, checkexit = TRUE)
}

# mothur make.shared
# here commands also differ because different clustering methods (above) produce different output files

if (datatype == "other") {
    makeshared.options <- paste("make.shared(list=file.opti_mcc.list, count=picked.count_table")
    makeshared.options <- paste(makeshared.options, ", label=", cutoff, ")", sep = "")
    documentCommand(makeshared.options)
    write(makeshared.options, "makeshared.mth", append = FALSE)
    command <- paste(binary, "makeshared.mth", "> log_makeshared.txt")
    system(command)
}

if (datatype == "its") {
    makeshared.options <- paste("make.shared(list=file.agc.list, count=picked.count_table") # produces file.agc.shared
    makeshared.options <- paste(makeshared.options, ", label=", cutoff, ")", sep = "")
    documentCommand(makeshared.options)
    write(makeshared.options, "makeshared.mth", append = FALSE)
    command <- paste(binary, "makeshared.mth", "> log_makeshared.txt")
    system(command)
}

# mothur classify.otu

if (datatype == "other") {
    classifyotu.options <- paste("classify.otu(list=file.opti_mcc.list, count=picked.count_table, taxonomy=sequences-taxonomy-assignment.txt")
    classifyotu.options <- paste(classifyotu.options, ", label=", cutoff, ")", sep = "")
    documentCommand(classifyotu.options)
    write(classifyotu.options, "classifyotu.mth", append = FALSE)
    command <- paste(binary, "classifyotu.mth", "> log_classifyotu.txt")
    system(command)
}

if (datatype == "its") {
    classifyotu.options <- paste("classify.otu(list=file.agc.list, count=picked.count_table, taxonomy=sequences-taxonomy-assignment.txt")
    classifyotu.options <- paste(classifyotu.options, ", label=", cutoff, ")", sep = "")
    documentCommand(classifyotu.options)
    write(classifyotu.options, "classifyotu.mth", append = FALSE)
    command <- paste(binary, "classifyotu.mth", "> log_classifyotu.txt")
    system(command)
}

# phenodata file
# based on mothur-classify-counttable.R

# read the data and tabulate it
pick <- read.table("picked.count_table", header = T, sep = "\t")
tax <- read.table("sequences-taxonomy-assignment.txt", header = F, sep = "\t")
dat <- merge(pick, tax, by.x = "Representative_Sequence", by.y = "V1")
dat$V2 <- gsub(".[[:digit:]]{1,}.?[[:digit:]]?)", "", as.character(dat$V2))

# cut taxonomic names
# based on default assumption of mothur-classify-counttable.R
# (i.e. that cutlevel = 0)
dat$newnames <- dat$V2

# set up the final result
data_start <- which(colnames(dat) == "total") + 1
data_end <- which(colnames(dat) == "V2") - 1
names_col <- which(colnames(dat) == "newnames")

# same manipulations here
dat <- dat[, c(names_col, data_start:data_end)]
datm <- melt(dat)
a <- aggregate(datm$value, list(datm$newnames, datm$variable), function(x) sum(x, na.rm = T))
b <- dcast(a, Group.2 ~ Group.1)
rownames(b) <- b$Group.2
b <- b[, -1]

tab <- b

# write phenodata.tsv
write.table(
    data.frame(
        sample = rownames(tab),
        chiptype = "NGS",
        group = rep("", length(rownames(tab)))
    ),
    "phenodata.tsv",
    col.names = T,
    row.names = F,
    sep = "\t",
    quote = F
)
