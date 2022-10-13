# TOOL cutadapt.R: "Remove primers and adpters with Cutadapt" (Given a tar package of fastq files, this tool removes the primer and adapter sequences.)
# INPUT OPTIONAL reads.tar: "Tar package containing the FASTQ files" TYPE GENERIC
# OUTPUT OPTIONAL adapters_removed.tar
# OUTPUT report.txt
# PARAMETER paired: "Is the data paired end or single end reads" TYPE [paired, single] DEFAULT single (Are all the reads paired end, so one forward and one reverse FASTQ file for one sample. If single end reads,use only those forward parameters.)
# PARAMETER OPTIONAL adapter5: "Forward or the 5' adapter to be trimmed" TYPE STRING
# PARAMETER OPTIONAL adapter3: "The 3' adapter to be trimmed" TYPE STRING
# RUNTIME R-4.1.1-asv

# ES 30.9.2022
# INPUT OPTIONAL file.fastq: "fastq file" TYPE GENERIC


source(file.path(chipster.common.path,"tool-utils.R"))
source(file.path(chipster.common.path,"zip-utils.R"))

# binary
binary <- c(file.path(chipster.tools.path,"python-3.8.11","bin","cutadapt"))
version <- system(paste(binary,"--version"),intern = TRUE)
documentVersion("Cutadapt",version)

if (adapter3=="" && adapter5==""){
  stop(paste('CHIPSTER-NOTE: ',"You need to give at least one adapter sequence to be removed"))
}

#check out if the file is compressed and if so unzip it
unzipIfGZipFile("reads.tar")

# Read the contents of the tar file into a list
system("tar tf reads.tar > tar.contents")
file.list <- scan("tar.contents",what = "",sep = "\n")

# Check that the input is a valid tar file
if (length(file.list) == 0) {
  stop(paste('CHIPSTER-NOTE: ',"It seems your input file is not a valid Tar package. Please check your input file."))
}
# Check if tar package contains folders
if (grepl("/",file.list[1])) {
  stop(paste('CHIPSTER-NOTE: ',"It seems your Tar package contains folders. The FASTQ files need to be in the root of the package, not in subfolders."))
}

# Make input and output folders 
system("mkdir input_folder")
system("mkdir output_folder")

# untar the tar package to input_folder
untar("reads.tar", exdir = "input_folder")

#list the full file names with the relative path
filenames <- sort(list.files("input_folder", full.names=TRUE))

# take out the names without the relative path
#short_names <- sub('\\..*', '', basename(filenames)) #take everything before first .
sample.names <- sapply(strsplit(basename(filenames), "_"), `[`, 1) #edit


#########################
# if paired end reads then need the forward and reverse reads in different places
if (paired=="paired"){

  # check if the lenght of files in the input_folder is even, else error
  number <- length(filenames)%%2
  if (number != 0){
      stop(paste('CHIPSTER-NOTE: ',"It seems that some of your fastq files doesn`t have a pair"))
      }
  # sorted filenames, allways forward,reverse,forward....
  forward <- seq(1,length(filenames),by=2)
  fnFs <- filenames[forward]
  fnFs.cut <- file.path("output_folder",paste0(basename(fnFs),"-cut-F.fastq.gz"))

  reverse <- seq(2,length(filenames), by=2)
  fnRs <- filenames[reverse]
  fnRs.cut<- file.path("output_folder",paste0(basename(fnRs),"-cut-R.fastq.gz"))
}else{
  #fnFs <- filenames
  cutreads <- file.path("output_folder",paste0(sample.names,"_cut.fastq.gz"))
  # Make output files to output_folder for every input file
}
######################## lets run cutadapt 

# use the parameters 3' and 5' and make the command, put the report to txt file, use --rc to check for reverse complements
# error message already sent if neither adapters selected, if paired use also R2.flag variable

if (adapter3==""){
  R1.flags <- paste("-g", adapter5)
  if (paired=="paired"){ #for reverse reads make a reverse complement
    R2.flags <- paste("-G", dada2:::rc(adapter5))
  }
}else if(adapter5==""){
  R1.flags <- paste("-a", adapter3)
  if (paired=="paired"){ #for reverse reads make a reverse complement
    R2.flags <- paste("-A", dada2:::rc(adapter3))}
}else{
  R1.flags <- paste("-g", adapter5, "-a", adapter3)   #dada2:::rc()
  if (paired=="paired"){ #for reverse reads make a reverse complement
    R2.flags <- paste("-G", dada2:::rc(adapter5), "-A", dada2:::rc(adapter3))
}}

#if paired use also R2.flags and 2 input files at once
if (paired =="single"){
  x<-1
  for (file in filenames){
    command <- paste(binary, R1.flags,"--rc","-n", 2,"-o",cutreads[x], file, "> report.txt")
    x=x+1
    system(command)
  }
}else{ #paired lets have fun
  x<-1
  for (file in fnFs){
    command <- paste(binary, R1.flags, R2.flags,"--rc","-n", 2,"-o",fnFs.cut[x],"-p", fnRs.cut[x], file, fnRs[x], "> report.txt")
    x=x+1
    system(command)
}
}



# make a tar package from the output folder
system("cd output_folder && tar cf ../adapters_removed.tar *")