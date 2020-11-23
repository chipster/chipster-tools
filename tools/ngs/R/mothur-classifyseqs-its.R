# TOOL mothur-classifyseqs-its.R: "Classify sequences ITS" (Explanation)
# INPUT a.fasta: "Aligned reads in FASTA file" TYPE FASTA
# INPUT a.count_table: "Count table" TYPE MOTHUR_COUNT
# INPUT OPTIONAL own_reference.fasta: "Reference FASTA" TYPE FASTA
# INPUT OPTIONAL own_reference.tax: "Reference taxonomy file" TYPE GENERIC
# OUTPUT OPTIONAL sequences-taxonomy-assignment.txt
# OUTPUT OPTIONAL classification-summary.tsv
# OUTPUT OPTIONAL picked.fasta.gz 
# OUTPUT OPTIONAL picked.count_table
# OUTPUT OPTIONAL picked-summary.tsv
# PARAMETER reference: "Reference" TYPE [own: "own reference", UNITEv8_sh_dynamic, UNITEv8_sh_99, UNITEv8_sh_97, UNITEv8_sh_dynamic_s, UNITEv8_sh_99_s, UNITEv8_sh_97_s] default own (Refernce to use.)
# PARAMETER OPTIONAL iters: "Number of iterations" TYPE INTEGER FROM 10 TO 1000 DEFAULT 100 (How many iterations to do when calculating the bootstrap confidence score for your taxonomy.)
# PARAMETER OPTIONAL remove.chloroplast: "Remove taxon Chloroplast" TYPE [yes, no] DEFAULT yes (Remove taxon Chloroplast.)
# PARAMETER OPTIONAL remove.mitochondria: "Remove taxon Mitochondria" TYPE [yes, no] DEFAULT yes (Remove taxon Mitochondria.)
# PARAMETER OPTIONAL remove.archaea: "Remove taxon Archaea" TYPE [yes, no] DEFAULT yes (Remove taxon Archaea.)
# PARAMETER OPTIONAL remove.eukaryota: "Remove taxon Eukaryota" TYPE [yes, no] DEFAULT yes (Remove taxon Eukaryota.)
# PARAMETER OPTIONAL remove.unknown: "Remove taxon unknown" TYPE [yes, no] DEFAULT yes (Remove taxon unknown.)
# PARAMETER OPTIONAL remove.other: "Remove other lineages" TYPE STRING DEFAULT empty (List of other lineages to remove. You must use dots \(\".\"\) instead of semicolons \(\";\"\). Use dash \(\"-\"\) to separate taxons. For example: Bacteria.Firmicutes.-Bacteria.Bacteroidetes. )

# OUTPUT OPTIONAL log.txt
# PARAMETER OPTIONAL reference: "Reference" TYPE [bacterial: "bacterial subset of Silva db", full: "whole Silva db"] DEFAULT bacterial (Silva reference set to use.)


# EK 18.06.2013
# JTT 28.8.2013 count table and phenodata added
# ML 21.12.2016 update (new Silva version)
# ML 14.3.2017 reference option (bacterial vs whole)
# ML 23.3.2017 detach the last steps to another tool (mothur-classify-counttable.R), add iters-parameter and remove.lineage option
# EK 22.8.2018 updated Silva to v 132, removed bacterial option, added processors-parameter, and made link in the working directory to the reference files
# EK 11.05.2020 Zip output fasta

source(file.path(chipster.common.path,"tool-utils.R"))
source(file.path(chipster.common.path,"zip-utils.R"))

# check out if the file is compressed and if so unzip it
unzipIfGZipFile("a.fasta")

# binary
binary <- c(file.path(chipster.tools.path,"mothur","mothur"))
version <- system(paste(binary,"--version"),intern = TRUE)
documentVersion("Mothur",version)

# if (reference=="bacterial"){
# bacterial references (Silva v102):
# 	data.path <- c(file.path(chipster.tools.path, "mothur-silva-reference", "silva.bacteria"))
# 	template.path <- c(file.path(data.path, "silva.bacteria.fasta"))
# 	taxonomy.path <- c(file.path(data.path, "silva.bacteria.silva.tax"))
# copy to working dir because mothur generates more files to the same directory
# 	system(paste("ln -s ", template.path, " silva.bacteria.fasta"))
# 	system(paste("ln -s ", taxonomy.path, " silva.bacteria.silva.tax"))
# 	template.path <- "silva.bacteria.fasta"
# 	taxonomy.path <- "silva.bacteria.silva.tax"
# }
if (reference == "own") {
  reference_path <- "own_reference.fasta"
  taxonomy.path <- "own_reference.tax"
} else {
  # Check data path when new tools-bin ready
  data.path <- c(file.path(chipster.tools.path,"mothur-its-reference","unite"))
  reference.file <- paste(reference,".fasta",sep = "")
  taxonomy.file <- paste(reference,".tax",sep = "")
  reference.path <- c(file.path(data.path,reference.file))
  taxonomy.path <- c(file.path(data.path,taxonomy.file))
  # copy to working dir because mothur generates more files to the same directory
  system(paste("ln -s ",reference.path,reference.file))
  system(paste("ln -s ",taxonomy.path,taxonomy.file))
  reference.path <- reference.file
  taxonomy.path <- taxonomy.file
}

# batch file
# write(paste("classify.seqs(fasta=a.fasta, iters=1000, template=", template.path, ", taxonomy=", taxonomy.path, ")", sep=""), "batch.mth", append=F)
classifyseqs.options <- paste("classify.seqs(fasta=a.fasta, count=a.count_table, ")
classifyseqs.options <- paste(classifyseqs.options,", iters=",iters,", refrence=",reference.path,", taxonomy=",taxonomy.path,", processors=",chipster.threads.max,")",sep = "")
documentCommand(classifyseqs.options)
write(classifyseqs.options,"batch.mth",append = FALSE)
# command
command <- paste(binary,"batch.mth","> log.txt 2>&1")

# run
system(command)

# Output File Names: 
# a.silva.wang.taxonomy
# a.silva.wang.tax.summary


## test
write("get.current()","batch2.mth",append = FALSE)
system(paste(binary,"batch2.mth",">> log.txt 2>&1"))

system("ls -l >> log.txt")
# Postprocess output
if (reference == "silva") {
  system("mv a.nr_v132.wang.taxonomy sequences-taxonomy-assignment.txt")
  system("mv a.nr_v132.wang.tax.summary classification-summary.tsv")
}
# if (reference=="bacterial"){
# 	system("mv a.silva.wang.taxonomy sequences-taxonomy-assignment.txt")
# 	system("mv a.silva.wang.tax.summary classification-summary.tsv")
# }

# batch file 2: remove lineage, if the taxons to remove were listed
toremove <- ""

if (remove.chloroplast == "yes") {
  toremove <- paste(toremove,"Chloroplast",sep = "-")
}
if (remove.mitochondria == "yes") {
  toremove <- paste(toremove,"Mitochondria",sep = "-")
}
if (remove.archaea == "yes") {
  toremove <- paste(toremove,"Archaea",sep = "-")
}
if (remove.eukaryota == "yes") {
  toremove <- paste(toremove,"Eukaryota",sep = "-")
}
if (remove.unknown == "yes") {
  toremove <- paste(toremove,"unknown",sep = "-")
}
if (remove.other != "empty") {
  # Change periods to semicolons. Semicolons are not accepted in STRING input.
  remove.other <- gsub(".",";",remove.other,fixed = TRUE)
  toremove <- paste(toremove,remove.other,sep = "-")
}

if (toremove != "") {
  # Remove leading dash
  toremove <- substring(toremove,2)
  removelineage.options <- paste("remove.lineage(fasta=a.fasta")
  if (file.exists("a.count_table")) {
    removelineage.options <- paste(removelineage.options,", count=a.count_table")
  }
  removelineage.options <- paste(removelineage.options,", taxonomy=sequences-taxonomy-assignment.txt, taxon=",toremove,")",sep = "")
  documentCommand(removelineage.options)
  write(removelineage.options,"batch.mth",append = FALSE)

  # command
  command <- paste(binary,"batch.mth",">> log.txt 2>&1")

  # run
  system(command)

  # Rename output files
  system("mv a.pick.fasta picked.fasta")
  if (file.exists("a.pick.count_table")) {
    system("mv a.pick.count_table picked.count_table")
  }

  # batch file 3: classify.seqs again

  # write(paste("classify.seqs(fasta=a.fasta, iters=1000, template=", template.path, ", taxonomy=", taxonomy.path, ")", sep=""), "batch.mth", append=F)
  classifyseqs.options <- paste("classify.seqs(fasta=picked.fasta")
  if (file.exists("picked.count_table")) {
    classifyseqs.options <- paste(classifyseqs.options,", count=picked.count_table")
  }
  classifyseqs.options <- paste(classifyseqs.options,", processors=",chipster.threads.max,", iters=",iters,", template=",template.path,", taxonomy=",taxonomy.path,")",sep = "")
  documentCommand(classifyseqs.options)
  write(classifyseqs.options,"batch.mth",append = FALSE)
  # command
  command <- paste(binary,"batch.mth",">> log.txt 2>&1")
  # run
  system(command)

  ## test
  write("get.current()","batch2.mth",append = FALSE)
  system(paste(binary,"batch2.mth",">> log.txt 2>&1"))


  # Postprocess output
  if (reference == "silva") {
    system("mv picked.nr_v132.wang.taxonomy sequences-taxonomy-assignment.txt")
    system("mv picked.nr_v132.wang.tax.summary classification-summary.tsv")
  }
  # if (reference=="bacterial"){
  # 	system("mv picked.silva.wang.taxonomy sequences-taxonomy-assignment.txt")
  # 	system("mv picked.silva.wang.tax.summary classification-summary.tsv")
  # }

  # batch file 3: summary
  summaryseqs.options <- paste("summary.seqs(fasta=picked.fasta, count=picked.count_table)")
  documentCommand(summaryseqs.options)
  write(summaryseqs.options,"summary.mth",append = FALSE)

  # command 3
  command3 <- paste(binary,"summary.mth","> log_raw.txt")

  # run
  system(command3)

  # zip output fasta
  system("gzip picked.fasta")

  # Post process output
  system("grep -A 10 Start log_raw.txt > picked-summary2.tsv")
  # Remove one tab to get the column naming look nice:
  system("sed 's/^		/	/' picked-summary2.tsv > picked-summary.tsv")
}
