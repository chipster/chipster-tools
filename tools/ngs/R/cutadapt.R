# TOOL cutadapt.R: "Remove primers and adapters with Cutadapt" (Given a tar package of FASTQ files, this tool tries to remove the primer and adapter sequences given in parameters tab. This tool is based on the tool Cutadapt )
# INPUT OPTIONAL reads.tar: "Tar package containing the FASTQ files" TYPE GENERIC
# INPUT OPTIONAL input_list.txt: "List of FASTQ files by sample" TYPE GENERIC (If the FASTQ files are not assigned into samples correctly, you can give a file containing this information. Check instructions from manual)
# OUTPUT adapters_removed.tar
# OUTPUT report.txt
# OUTPUT OPTIONAL samples.fastqs.txt
# PARAMETER paired: "Is the data paired end or single end reads" TYPE [paired, single] DEFAULT paired (If your reads are paired end, then the adapters from the reverse files will be removed by removing the reverse complement of the 3' and 5' adapters given as parameters.)
# PARAMETER OPTIONAL adapter5: "The 5' adapter:" TYPE STRING (Give here the 5 end adapter/primer)
# PARAMETER OPTIONAL adapter3: "The 3' adapter:" TYPE STRING (Give here the 3 end adapter/primer)
# RUNTIME R-4.1.1-asv

# ES 30.9.2022
# ES 28.12.2022 added samples.fastqs.txt output and input file
# added multi-core support 20.10.2022
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

# without path to make the samples fastqs txt file
txt_filenames <- sort(list.files("input_folder", full.names=FALSE))

# take out the names without the relative path
#short_names <- sub('\\..*', '', basename(filenames)) #take everything before first .
sample.names <- sapply(strsplit(basename(filenames), "_"), `[`, 1) 


#########################
# if paired end reads, put the forward and reverse files to separate folders
if (paired=="paired"){

# check if the lenght of files in the input_folder is even, else error
  number <- length(filenames)%%2
  if (number != 0){
      stop(paste('CHIPSTER-NOTE: ',"It seems that some of your fastq files doesn`t have a pair"))
      }

  # Use input_list if provided. Else use the sorted filenames from the tar package
  if (fileOk("input_list.txt")) {
    txt_filenames <- c() # name of the files in the txt file
    sample.names <- c() #sample names
    input <- readLines("input_list.txt") #read the file

    # take out the filenames and sample name and put them to one vector, those are separeted with '\t'
    for (row in input){
      sample <- strsplit(row,'\t',fixed=TRUE)
      sample_name<- trimws(sample[[1]][1])
      sample.names <- c(sample.names, sample_name) # sample names first row
      txt_filenames <- c(txt_filenames, trimws(sample[[1]][2])) # forward read
      txt_filenames <- c(txt_filenames, trimws(sample[[1]][3])) # reverse read
    }
    # if everything fine change the filenames variable and use it, add also the folder name: input_folder/
    if (length(txt_filenames) != length(filenames)){
      print(txt_filenames, filenames)
      stop(paste('CHIPSTER-NOTE: ',"It seems that the list of FASTQ files .txt file has different amount of filenames than the .tar package"))
    }else{
   
      filenames <- paste0("input_folder/",txt_filenames) # now the input files in correct order
    }
  }

  # put filenames with full paths to separate folders
  # forward files:
  forward <- seq(1,length(filenames),by=2)
  fnFs <- filenames[forward]

  # sample names unique 
  unique_samples <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
  fnFs.cut <- file.path("output_folder",paste0(unique_samples,"_cut_F.fastq.gz"))

  # reverse files:
  reverse <- seq(2,length(filenames), by=2)
  fnRs <- filenames[reverse]
  fnRs.cut<- file.path("output_folder",paste0(unique_samples, "_cut_R.fastq.gz"))

}else{ # single end reads
  #fnFs <- filenames
  cutreads <- file.path("output_folder",paste0(sample.names,"_cut.fastq.gz"))
  # Make output files to output_folder for every input file
}
######################## lets run cutadapt 

# use the parameters 3' and 5' and make the command, put the report to txt file, use --rc to check for reverse complements
# error message already sent if neither adapters selected, if paired use also R2.flag variable
# create a report file to combine all reports 
file.create("report.txt")

if (adapter3==""){
  R1.flags <- paste("-g", adapter5)
  if (paired=="paired"){ #for reverse reads make a reverse complement
    R2.flags <- paste("-A", dada2:::rc(adapter5))
  }
}else if(adapter5==""){
  R1.flags <- paste("-a", adapter3)
  if (paired=="paired"){ #for reverse reads make a reverse complement
    R2.flags <- paste("-G", dada2:::rc(adapter3))} #s
}else{
  R1.flags <- paste("-g", adapter5, "-a", adapter3)   #dada2:::rc()
  if (paired=="paired"){ #for reverse reads make a reverse complement and other way around
    #R2.flags <- paste("-G", dada2:::rc(adapter5), "-A", dada2:::rc(adapter3))
    R2.flags <- paste("-G", dada2:::rc(adapter3), "-A", dada2:::rc(adapter5))
}}

#sink(file="report.txt")
#if paired use also R2.flags and 2 input files at once
if (paired =="single"){
  x<-1
  for (file in filenames){
    command <- paste(binary, R1.flags,"--rc","-n", 2,"-j", as.integer(chipster.threads.max),"-o",cutreads[x], file, "> report2.txt")
    x=x+1
    system(command)
    system("cat report2.txt >> report.txt")
  }
}else{ #paired 
  x<-1
  for (file in fnFs){
    command <- paste(binary, R1.flags, R2.flags,"-n", 2,"-j", as.integer(chipster.threads.max),"-o",fnFs.cut[x],"-p", fnRs.cut[x], file, fnRs[x] , "> report2.txt")
    x=x+1
    system(command)
    #rows <- readLines("report2.txt")
    #for (row in rows) {
     # write(row, file = "report.txt", append=TRUE)
   #}
    #write("\n", file = "report.txt", append=TRUE)
    system("cat report2.txt >> report.txt")
    
}

# make samples.fastqs.txt file to show how FASTQ files were assigned to samples, every second file forward if paired end reads
file.create("samples.fastqs.txt")
x <-1
for (name in unique_samples){
  line <- paste(name,'\t',txt_filenames[x],'\t',txt_filenames[x+1])
  x=x+2
  write(line, file="samples.fastqs.txt", append=TRUE)
}
}
#sink()

# make a tar package from the output folder
system("cd output_folder && tar cf ../adapters_removed.tar *")