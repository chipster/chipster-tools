# TOOL mothur-phyloseq-inputs.R: "Generate input files for phyloseq" (This tool generates input files for the R package phyloseq, using the Mothur tools dist.seqs, cluster, make.shared and classify.otu. As part of input file generation, sequences are clustered into operational taxonomic units \(OTU\) based on a user-defined sequence dissimilarity threshold. Note that sample groups in the phenodata file should be indicated using characters \(numeric labels are currently not supported by downstream analyses\).)
# INPUT file.fasta: "FASTA file" TYPE GENERIC
# INPUT picked.count_table: "Mothur count file" TYPE MOTHUR_COUNT
# INPUT sequences-taxonomy-assignment.txt: "Sequences taxonomy assignment file" TYPE GENERIC
# OUTPUT META phenodata.tsv
# OUTPUT file.opti_mcc.shared
# OUTPUT log_distseqs.txt
# OUTPUT OPTIONAL file.opti_mcc.0.05.cons.taxonomy
# OUTPUT OPTIONAL file.opti_mcc.0.04.cons.taxonomy
# OUTPUT OPTIONAL file.opti_mcc.0.03.cons.taxonomy
# OUTPUT OPTIONAL file.opti_mcc.0.02.cons.taxonomy
# OUTPUT OPTIONAL file.opti_mcc.0.01.cons.taxonomy
# PARAMETER cutoff: "Cutoff" TYPE [0.05, 0.04, 0.03, 0.02, 0.01] DEFAULT 0.03 (Dissimilarity threshold for OTU clustering, e.g. a cut-off value of 0.03 corresponds to 97% similarity)

# reshape2 library
library(reshape2)

# check out if the file is compressed and if so unzip it
source(file.path(chipster.common.path,"tool-utils.R"))
source(file.path(chipster.common.path,"zip-utils.R"))
unzipIfGZipFile("file.fasta")

# binary
binary <- c(file.path(chipster.tools.path,"mothur","mothur"))
version <- system(paste(binary,"--version"),intern = TRUE)
documentVersion("Mothur",version)

# mothur dist.seqs
# produces:
# file.dist

distseqs.options <- paste("dist.seqs(fasta=file.fasta")
distseqs.options <- paste(distseqs.options,", processors=",chipster.threads.max,sep = "")
distseqs.options <- paste(distseqs.options,", cutoff=",cutoff,")",sep = "")
documentCommand(distseqs.options)
write(distseqs.options,"distseqs.mth",append = FALSE)
command <- paste(binary,"distseqs.mth","> log_distseqs.txt")
# system(command)
runExternal(command, checkexit = TRUE)

# mothur cluster
# produces:
# file.opti_mcc.list
# file.opti_mcc.steps
# file.opti_mcc.sensspec

cluster.options <- paste("cluster(column=file.dist, count=picked.count_table")
cluster.options <- paste(cluster.options,", cutoff=",cutoff,")",sep = "")
documentCommand(cluster.options)
write(cluster.options,"cluster.mth",append = FALSE)
command <- paste(binary,"cluster.mth","> log_cluster.txt")
# system(command)
runExternal(command, checkexit = TRUE)

# mothur make.shared
# produces:
# file.opti_mcc.shared

makeshared.options <- paste("make.shared(list=file.opti_mcc.list, count=picked.count_table")
makeshared.options <- paste(makeshared.options,", label=",cutoff,")",sep = "")
documentCommand(makeshared.options)
write(makeshared.options,"makeshared.mth",append = FALSE)
command <- paste(binary,"makeshared.mth","> log_makeshared.txt")
system(command)

# mothur classify.otu
# produces:
# file.opti_mcc.[cutoff].cons.taxonomy
# file.opti_mcc.[cutoff].cons.tax.summary

classifyotu.options <- paste("classify.otu(list=file.opti_mcc.list, count=picked.count_table, taxonomy=sequences-taxonomy-assignment.txt")
classifyotu.options <- paste(classifyotu.options,", label=",cutoff,")",sep = "")
documentCommand(classifyotu.options)
write(classifyotu.options,"classifyotu.mth",append = FALSE)
command <- paste(binary,"classifyotu.mth","> log_classifyotu.txt")
system(command)

# phenodata file
# based on mothur-classify-counttable.R

# read the data and tabulate it
pick <- read.table("picked.count_table",header = T,sep = "\t")
tax <- read.table("sequences-taxonomy-assignment.txt",header = F,sep = "\t")
dat <- merge(pick,tax,by.x = "Representative_Sequence",by.y = "V1")
dat$V2 <- gsub(".[[:digit:]]{1,}.?[[:digit:]]?)","",as.character(dat$V2))

# cut taxonomic names
# based on default assumption of mothur-classify-counttable.R
# (i.e. that cutlevel = 0)
dat$newnames <- dat$V2

# set up the final result
data_start <- which(colnames(dat) == "total") + 1
data_end <- which(colnames(dat) == "V2") - 1
names_col <- which(colnames(dat) == "newnames")

# same manipulations here
dat <- dat[,c(names_col,data_start:data_end)]
datm <- melt(dat)
a <- aggregate(datm$value,list(datm$newnames,datm$variable),function(x) sum(x,na.rm = T))
b <- dcast(a,Group.2 ~ Group.1)
rownames(b) <- b$Group.2
b <- b[,-1]

tab <- b

# write phenodata.tsv
write.table(data.frame(sample = rownames(tab),
    chiptype = "NGS",
    group = rep("",length(rownames(tab)))),
    "phenodata.tsv",
    col.names = T,
    row.names = F,
    sep = "\t",
    quote = F)
