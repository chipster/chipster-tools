# TOOL fastx-trimmer.R: "Trim reads with FastX" (Trims reads to a user-specified length. Input can be one or several FASTQ files, and the files can be in a tar package which can be gzipped. If the input files are in a tar package, then also the result files are put in a tar package. This tool is based on the FASTA/Q Trimmer tool of the FASTX package.)
# INPUT input.file TYPE GENERIC
# OUTPUT OPTIONAL trimmed.fq.gz
# OUTPUT OPTIONAL trimmed.tar
# PARAMETER first: "First base to keep" TYPE INTEGER FROM 1 TO 100 DEFAULT 1 (First base to keep.)
# PARAMETER last: "Last base to keep" TYPE INTEGER FROM 1 TO 300 DEFAULT 75 (Last base to keep.)
# PARAMETER quality.format: "Quality value format used" TYPE [sanger: Sanger, illuminaold: "Illumina GA v1.3-1.5"] DEFAULT sanger (What quality encoding is used in your FASTQ file. Select Sanger if your data comes from Illumina 1.8 or later, SOLiD or 454.)




# EK 17.6.2011
# AMS 11.3.2014, gzip fastq outputs

source(file.path(chipster.common.path,"tool-utils.R"))
source(file.path(chipster.common.path, "zip-utils.R"))



# Binary
binary <- c(file.path(chipster.tools.path, "fastx", "bin", "fastx_trimmer"))
version <- system(paste(binary,"-h | sed -n 2p"),intern = TRUE)
documentVersion("fastx_trimmer",version)

# Check if input is a tar file
isTar <- grepl("POSIX tar", system("file input.file", intern=TRUE))

if (isTar){
  # Input is a tar file
  system("mkdir input_folder")
  system("mkdir output_folder")
  system("cd input_folder && tar xf ../input.file")
  system("cd input_folder && gunzip *.gz")
  filenames <- list.files("input_folder")
  for (f in filenames){
    input_fastq <- paste("input_folder/",f, sep="")
    output_base <- strip_name(f)
    output_fastq <- paste("output_folder/",output_base, "_trimmed.fq", sep="")

    # Command
    quality.scale <- ifelse(quality.format == "sanger", "-Q 33", "")
    command <- paste(binary, "-f", first, "-l", last, quality.scale, "-i", input_fastq, "-o", output_fastq)
    #documentCommand(command)
    system(command)
  }
  # gzip all output FASTQ files
  system("gzip output_folder/*.fq")
  # Make a tar package.
  system("cd output_folder && tar cf ../trimmed.tar *")

}else{
  # Input is FASTQ file

  # Check out if the file is compressed and if so unzip it
  unzipIfGZipFile("input.file")

  # Command
  quality.scale <- ifelse(quality.format == "sanger", "-Q 33", "")
  command <- paste(binary, "-f", first, "-l", last, quality.scale, "-i input.file -o trimmed.fq")
  documentCommand(command)

  # Run
  system(command)
  system("gzip *.fq")

  # Determine base name
  inputnames <- read_input_definitions()
  basename <- strip_name(inputnames$input.file)

  # Make a matrix of output names
  outputnames <- matrix(NA,nrow = 1,ncol = 2)
  outputnames[1,] <- c("trimmed.fq.gz",paste(basename,"_trimmed.fq.gz",sep = ""))

  # Write output definitions file
  write_output_definitions(outputnames)
}
