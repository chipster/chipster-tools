# TOOL test-mothur.R: "Test-Mothur"
# INPUT file.fasta: "FASTA file" TYPE GENERIC
# INPUT final.count_table: "Mothur count file" TYPE MOTHUR_COUNT
# INPUT sequences-taxonomy-assignment.txt: "Sequences taxonomy assignment file" TYPE GENERIC
# OUTPUT OPTIONAL final.unique.list
# OUTPUT OPTIONAL final.asv.shared
# OUTPUT OPTIONAL final.asv.list
# OUTPUT META phenodata.tsv
# OUTPUT OPTIONAL final.unique.shared
# OUTPUT OPTIONAL log_cluster.txt
# OUTPUT OPTIONAL log_distseqs.txt
# OUTPUT OPTIONAL log_makeshared.txt
# OUTPUT OPTIONAL log_classifyotu.txt
# OUTPUT OPTIONAL final.unique.0.03.cons.taxonomy
# OUTPUT OPTIONAL final.asv.asv.cons.taxonomy
# OUTPUT OPTIONAL final.asv.asv.cons.tax.summary


# check out if the file is compressed and if so unzip it
source(file.path(chipster.common.path,"tool-utils.R"))
source(file.path(chipster.common.path,"zip-utils.R"))
unzipIfGZipFile("file.fasta")

# binary
binary <- c(file.path(chipster.tools.path,"mothur","mothur"))
version <- system(paste(binary,"--version"),intern = TRUE)
documentVersion("Mothur",version)

library(reshape2)

#distseqs.options <- paste("dist.seqs(fasta=file.fasta)") # dist.seqs produces file.dist
#distseqs.options <- paste(distseqs.options,", processors=",chipster.threads.max,sep = "")
#distseqs.options <- paste(distseqs.options,", cutoff=",cutoff,")",sep = "")
#documentCommand(distseqs.options)
#write(distseqs.options,"distseqs.mth",append = FALSE)
#command <- paste(binary,"distseqs.mth","> log_distseqs.txt")
#system(command)
#runExternal(command, checkexit = TRUE)



cluster.options <- paste("cluster(fasta=file.fasta, count=final.count_table, method=unique)") #column=file.dist
documentCommand(cluster.options)
write(cluster.options,"cluster.mth",append = FALSE)
command <- paste(binary,"cluster.mth","> log_cluster.txt")
system(command)
runExternal(command)

#makeshared.options <- paste("make.shared(list=final.unique.list, count=final.count_table)")
#makeshared.options <- paste(makeshared.options,", label=asv)",sep = "")
makeshared.options <- paste("make.shared(count=final.count_table, label=asv)")
documentCommand(makeshared.options)
write(makeshared.options,"makeshared.mth",append = FALSE)
command <- paste(binary,"makeshared.mth","> log_makeshared.txt")
system(command)

classifyotu.options <- paste("classify.otu(list=final.asv.list, count=final.count_table, taxonomy=sequences-taxonomy-assignment.txt, label=asv)")
documentCommand(classifyotu.options)
write(classifyotu.options,"classifyotu.mth",append = FALSE)
command <- paste(binary,"classifyotu.mth","> log_classifyotu.txt")
system(command)

# read the data and tabulate it
pick <- read.table("final.count_table",header = T,sep = "\t")
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
