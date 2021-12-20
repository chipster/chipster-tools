# TOOL subset-fastq.R: "Make a subset of FASTQ" (Returns a subsample of N reads from a FASTQ file. The input can be a single FASTQ file or a tar package containing multiple FASTQ files. Tool is based on the seqtk package. When using paired-end data, use same random seed to keep pairing.)
# INPUT input.file: "FASTQ file" TYPE GENERIC
# OUTPUT OPTIONAL subset.fastq.gz
# OUTPUT OPTIONAL subset.tar 
# PARAMETER n.seq: "Size of subset" TYPE INTEGER DEFAULT 100000 (Number of reads to return from the FASTQ file.)
# PARAMETER seed: "Random seed" TYPE INTEGER DEFAULT 11 (Random seed for the sampling. When using paired-end data, use same random seed to keep pairing.)

# AMS 14.5.2012
# EK 15.5.2012 added unzipping
# AMS 24.9.2014 added zipping the result file
# AMS 13.10.2014 changed to use seqtk


source(file.path(chipster.common.path, "zip-utils.R"))
source(file.path(chipster.common.path, "zip-utils.R"))

# binary
seqtk.binary <- file.path(chipster.tools.path, "seqtk", "seqtk")

seed.option <- paste ("-s",as.character(seed), sep="")

# Check if input is a tar file
isTar <- grepl("POSIX tar", system("file input.file", intern=TRUE))

if (isTar){
  # Input is a tar file
  system("mkdir input_folder")
  system("mkdir output_folder")
  system("cd input_folder && tar xf ../input.file")
  system("cd input_folder && gunzip *.gz")
  system("cd input_folder && ls -l")
  filenames <- list.files("input_folder")
  for (f in filenames){
    input_fastq <- paste("input_folder/",f, sep="")
    output_base <- strip_name(f)
    output_fastq <- paste("output_folder/",output_base, "_subset.fq", sep="")

    # Command
    command <- paste(seqtk.binary, "sample", seed.option, input_fastq, n.seq, ">", output_fastq)
    documentCommand(command)
    runExternal(command)
  }  
  # gzip all output FASTQ files
  system("gzip output_folder/*.fq")
  # Make a tar package.
  system("cd output_folder && tar cf ../subset.tar *")

  # read input names
  inputnames <- read_input_definitions()
  base <- strip_name(inputnames$input.file)

  # Make a matrix of output names
  outputnames <- matrix(NA, nrow=1, ncol=2)
  outputnames[1,] <- c("subset.tar", paste(base, "_subset.tar", sep =""))

  # Write output definitions file
  write_output_definitions(outputnames)

}else{
  # command
  unzipIfGZipFile("reads.fastq")
  command <- paste(seqtk.binary, "sample", seed.option, "input.file", n.seq, "> subset.fastq")
  documentCommand(command)
  runExternal(command)
  system("gzip subset.fastq")

  # read input names
  inputnames <- read_input_definitions()
  base <- strip_name(inputnames$input.file)

  # Make a matrix of output names
  outputnames <- matrix(NA, nrow=1, ncol=2)
  outputnames[1,] <- c("subset.fastq.gz", paste(base, "_subset.fq.gz", sep =""))

  # Write output definitions file
  write_output_definitions(outputnames)
}