# TOOL mothur-alignseqs.R: "Align sequences to reference" (Given a fasta file of 16S rRNA sequences, this tool aligns them to the Silva reference alignment available on the server, or to a reference alignment supplied by the user. Please indicate the region of the reference alignment that your amplified region corresponds to. If you opt to use your own reference alignment, please make sure the input fasta files are correctly assigned in the parameters section! The speed of this tool depends on several things, please see the manual. This tool is based on the Mothur tool align.seqs.)
# INPUT reads.fasta: "FASTA file" TYPE FASTA
# INPUT OPTIONAL reference.fasta: "Custom reference FASTA file" TYPE FASTA
# INPUT OPTIONAL a.count_table: "Count table" TYPE MOTHUR_COUNT
# OUTPUT aligned.fasta.gz
# OUTPUT aligned-summary.tsv
# OUTPUT custom.reference.summary.tsv
# OUTPUT log.txt
# PARAMETER OPTIONAL reference: "Reference" TYPE [silva: "silva.nr_v138.1", own: "own reference in fasta format"] DEFAULT silva (Reference sequence alignment to use.)
# PARAMETER OPTIONAL start: "Start" TYPE INTEGER (Start point of your region of interest)
# PARAMETER OPTIONAL end: "End" TYPE INTEGER (End point of your region of interest)
# SLOTS 5


# EK 05.06.2013
# ML 21.12.2016 update (new Silva version)
# ML 4.1.2017 new, whole Silva reference
# ML 14.3.2017 reference option (bacterial vs whole)
# ML 15.3.2017 add pcr.seqs options
# EK 22.8.2018 updated Silva to v132, added processors parameter to pcr.seqs and align.seqs
# EK 19.4.2021 updated to Silva v138.1

# PARAMETER OPTIONAL keepdots: "Remove leading and trailing dots" TYPE [yes, no] DEFAULT yes (Remove leading and trailing dots.)


source(file.path(chipster.common.path,"tool-utils.R"))
source(file.path(chipster.common.path,"zip-utils.R"))

# check out if the file is compressed and if so unzip it
unzipIfGZipFile("reads.fasta")

# binary
binary <- c(file.path(chipster.tools.path,"mothur","mothur"))
version <- system(paste(binary,"--version"),intern = TRUE)
documentVersion("Mothur",version)


#if (reference=="bacterial"){
# new bacterial references:
#	data.path <- c(file.path(chipster.tools.path, "mothur-silva-reference", "silva.bacteria"))
#	template.path <- c(file.path(data.path, "silva.bacteria.fasta"))
#}
#if (reference=="full"){

# new whole references:
# data.path <- c(file.path(chipster.tools.path,"mothur-silva-reference", "mothur-silva-reference-whole"))
data.path <- c(file.path(chipster.tools.path,"mothur-silva-reference","silva"))
# template.path <- c(file.path(data.path,"silva.nr_v132.align"))
template.path <- c(file.path(data.path,"silva.nr_v138_1.align"))
#}

# create a symlink, because otherwise the modified reference will go to the reference folder
system(paste("ln -s ",template.path," template.fasta",sep = ""))


# batch file 1 -pcr.seqs -only if user determined end or start 

if (!is.na(start) | !is.na(end)) {

  pcrseqs.options <- ""
  if (reference == "own") {
    if (file.exists("reference.fasta")) {
      pcrseqs.options <- paste(pcrseqs.options,"pcr.seqs(fasta=reference.fasta, processors=",chipster.threads.max,", keepdots=F",sep = "")
    } else {
      stop('CHIPSTER-NOTE: If you choose to use your own reference, you need to give the fasta file for that as input!')
    }
  } else {
    # if using silva reference on the server:
    #pcrseqs.options <- paste(pcrseqs.options, "pcr.seqs(fasta=", template.path, sep="")	
    pcrseqs.options <- paste(pcrseqs.options,"pcr.seqs(fasta=template.fasta, processors=",chipster.threads.max,", keepdots=F",sep = "")
  }

  if (!is.na(start)) {
    pcrseqs.options <- paste(pcrseqs.options,", start=",start,sep = "")
  }
  if (!is.na(end)) {
    pcrseqs.options <- paste(pcrseqs.options,", end=",end,sep = "")
  }
  #	if (keepdots=="yes"){
  #		pcrseqs.options <- paste(pcrseqs.options, ", keepdots=F", sep="")
  #	}
  #	if (keepdots=="no"){
  #		pcrseqs.options <- paste(pcrseqs.options, ", keepdots=T", sep="")
  #	}
  pcrseqs.options <- paste(pcrseqs.options,", processors=",chipster.threads.max,sep = "")
  pcrseqs.options <- paste(pcrseqs.options,")",sep = "")

  # Write batch file
  documentCommand(pcrseqs.options)
  write(pcrseqs.options,"pcrseq.mth",append = FALSE)
  # command
  command <- paste(binary,"pcrseq.mth","> log.txt 2>&1")
  # run
  system(command)

}

# rename the reference file as custom.reference.fasta 
if (file.exists("template.pcr.fasta")) {
  system("mv template.pcr.fasta custom.reference.fasta")
} else if (file.exists("reference.pcr.fasta")) {
  system("mv reference.pcr.fasta custom.reference.fasta")
} else if (file.exists("template.fasta")) {
  system("mv template.fasta custom.reference.fasta")
}



#  summary file from this step:
if (file.exists("custom.reference.fasta")) {
  # batch file 2
  summaryseqs.options <- paste("summary.seqs(fasta=custom.reference.fasta")
  summaryseqs.options <- paste(summaryseqs.options,", processors=",chipster.threads.max,")",sep = "")
  documentCommand(summaryseqs.options)
  write(summaryseqs.options,"summary.mth",append = FALSE)
  # command
  command2 <- paste(binary,"summary.mth","> log_raw.txt")
  # run
  system(command2)
  # Post process output
  system("grep -A 10 Start log_raw.txt > custom.reference.summary2.tsv")
  # Remove one tab to get the column naming look nice:
  system("sed 's/^		/	/' custom.reference.summary2.tsv > custom.reference.summary.tsv")
}


# batch file 2 -align.seqs
alignseqs.options <- paste("align.seqs(fasta=reads.fasta, processors=",chipster.threads.max,", reference=custom.reference.fasta)",sep = "")
documentCommand(alignseqs.options)
write(alignseqs.options,"batch.mth",append = FALSE)
# write(paste("align.seqs(fasta=reads.fasta, template=", template.path, ")", sep=""), "batch.mth", append=F)

# command
command <- paste(binary,"batch.mth",">> log.txt 2>&1")

# run
system(command)

system("mv reads.align aligned.fasta")

# batch file 3 -summary.seqs
summaryseqs.options <- paste("summary.seqs(fasta=aligned.fasta")
if (file.exists("a.count_table")) {
  summaryseqs.options <- paste(summaryseqs.options,", count=a.count_table")
}
summaryseqs.options <- paste(summaryseqs.options,", processors=",chipster.threads.max,")",sep = "")
documentCommand(summaryseqs.options)
write(summaryseqs.options,"summary.mth",append = FALSE)
# command
command2 <- paste(binary,"summary.mth","> log_raw.txt")

# run
system(command2)

# zip fasta
system("gzip aligned.fasta")

# Post process output
system("grep -A 10 Start log_raw.txt > aligned-summary2.tsv")
# Remove one tab to get the column naming look nice:
system("sed 's/^		/	/' aligned-summary2.tsv > aligned-summary.tsv")
