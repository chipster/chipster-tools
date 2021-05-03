# TOOL mothur-makegroup.R: "Combine FASTQ or FASTA files and make a group file" (Combines FASTQ or FASTA files of all samples to one file and creates a Mothur group file for it. Input can be one or more FASTQ or FASTA files, or a tar file. Input files can be gzipped. This tool is based on the Mothur tools make.group, merge.files and summary.seqs.)
# INPUT reads{...}: "read files" TYPE GENERIC
# OUTPUT sequences.fasta.gz
# OUTPUT sequences.groups
# OUTPUT sequences-summary.tsv

# EK 6.4.2021
# OUTPUT log1.txt


source(file.path(chipster.common.path, "tool-utils.R"))
source(file.path(chipster.common.path, "zip-utils.R"))

# check out if the file is compressed and if so unzip it
input.names <- read.table("chipster-inputs.tsv", header=F, sep="\t")

# Binary
fastx.binary <- c(file.path(chipster.tools.path, "fastx", "bin", "fastq_to_fasta"))

# Go through the input files. Uncompress if necessary. Untar if necessary. 
fasta_names <- ""
group_names <- ""
for (i in 1:nrow(input.names)) {
	unzipIfGZipFile(input.names[i,1])
  # Check if input is tar package
  isTar <- grepl("POSIX tar", system(paste("file", input.names[i,1]), intern=TRUE))
  # Input is a tar file
  if (isTar){
    tarlist <- untar(paste(input.names[i,1]), list = TRUE)
    untar(paste(input.names[i,1]))
    # go throug list
    for (j in 1:length(tarlist)){
      unzipIfGZipFile(tarlist[j])
      basename <- paste(strip_name(paste(tarlist[j])))
      fastaname <- paste(basename, ".fasta", sep = "")
      # If input is FASTQ file, convert to FASTA
      if (isFastq(paste(tarlist[j]))){
        command <- paste(fastx.binary, "-n -i", tarlist[j],"-o", fastaname)
        runExternal(command)        
      # If input is FASTA file, just make sure file name ends with .fasta  
      }else if(isFasta(paste(tarlist[j]))){
        file.rename(paste(tarlist[j]), fastaname)
      } else{
        stop(paste('CHIPSTER-NOTE: Input not a FASTA or FASTQ sequence file'))
      }
      fasta_names <- paste(fasta_names, fastaname, sep="-")
      group_names <- paste(group_names, basename, sep="-")
    }
  }
  # Input file is a sequense file
  else{
    basename <- paste(strip_name(paste(input.names[i,1])))
    fastaname <- paste(basename, ".fasta", sep = "")
    groupname <- paste(strip_name(paste(input.names[i,2])))
    # Input is a FASTQ file, convert to FASTA    
    if (isFastq(paste(input.names[i,1]))){      
      command <- paste(fastx.binary, "-n -i", input.names[i,1],"-o", fastaname)
      runExternal(command)
    # Input is FASTA file, just make sure file name ends with .fasta
    }else if(isFasta(paste(input.names[i,1]))){
      file.rename(paste(input.names[i,1]), fastaname)    
    }else{
      stop(paste('CHIPSTER-NOTE: Input not a FASTA or FASTQ sequence file'))
    }
    fasta_names <- paste(fasta_names, fastaname, sep="-")
    group_names <- paste(group_names, groupname, sep="-")
  }
}
# List files now have extra "-" at beginning that has to be removed
fasta_names <- sub(".","", fasta_names)
group_names <- sub(".","",group_names)

# Binary
binary <- c(file.path(chipster.tools.path,"mothur","mothur"))
version <- system(paste(binary,"--version"),intern = TRUE)
documentVersion("Mothur",version)

makegroup.options <- paste("make.group(fasta=",fasta_names, ", groups=",group_names, ")", sep = "")
documentCommand(makegroup.options)
write(makegroup.options,"makegroup.mth",append = FALSE)
command1 <- paste(binary,"makegroup.mth","> log1.txt")
system(command1)

# Rename groups file
filename <- list.files(path =".", pattern="*groups")
system(paste("mv", filename[1], "sequences.groups"))

# Merge fasta files
mergefiles.options <- paste("merge.files(input=", fasta_names, ", output=sequences.fasta)",  sep = "")
documentCommand(mergefiles.options)
write(mergefiles.options,"mergefiles.mth",append = FALSE)
command2 <- paste(binary,"mergefiles.mth","> log2.txt")
system(command2)

# Create summary file
summaryseqs.options <- paste("summary.seqs(fasta=sequences.fasta")
summaryseqs.options <- paste(summaryseqs.options,", processors=",chipster.threads.max,")",sep = "")
documentCommand(summaryseqs.options)
write(summaryseqs.options,"summary.mth",append = FALSE)
command <- paste(binary,"summary.mth","> log_raw.txt")
system(command)
# Post process output
system("grep -A 10 Start log_raw.txt > summary.tsv")
# Remove one tab to get the column naming look nice:
system("sed 's/^		/	/' summary.tsv > sequences-summary.tsv")

# Gzip output fasta
system("gzip sequences.fasta")

#Sys.sleep(120)
