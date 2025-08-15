# TOOL basespace.R: "Retrieve data from Illumina BaseSpace" (Retrieve data from Illumina BaseSpace. This tool requires that you have an access token for the bs client program. Please see the manual for how to obtain it.)
# OUTPUT OPTIONAL bs_data.tsv
# OUTPUT OPTIONAL bs_data.txt
# OUTPUT OPTIONAL bs_data.tar
# PARAMETER OPTIONAL name: "Name of dataset or project" TYPE STRING (Give the name of the dataset or project. This parameter is not needed if you just want to list your datasets or projects in Illumina BaseSpace.)
# PARAMETER action: "Action" TYPE [list_datasets: "List datasets", list_projects: "List projects", download_dataset: "Download dataset", download_project: "Download project", dir: "Display content of a dataset", info: "Display detailed information about a dataset"  ] DEFAULT list_datasets (Action to be performed.)
# PARAMETER apiserver: "API server" TYPE [api.basespace.illumina.com: "api.basespace.illumina.com"] DEFAULT api.basespace.illumina.com (Define the BaseSpace server to be used.)
# PARAMETER token: "Access token" TYPE STRING (Your personal Illumina BaseSpace access token.)

# OUTPUT OPTIONAL {...}.fastqc.gz: "FASTQ files"

# KM 25.02.2018

source(file.path(chipster.common.lib.path, "tool-utils.R"))

# Make a matrix of output names
outputnames <- matrix(NA, nrow = 1, ncol = 2)


bs.binary <- file.path(chipster.tools.path, "basespace/bin/bs")
# turn of cacheing
bs.command_start <- paste(bs.binary, " --api-server=https://", apiserver, " --access-token=", token, sep = "")

if (action == "list_datasets") {
  command.full <- paste(bs.command_start, 'list dataset -f csv | tr "," "\t" 1>>bs_data.tsv')
  runExternal(command.full)
  outputnames[1, ] <- c("bs_data.tsv", "basespace-datasets.tsv")
}

if (action == "list_projects") {
  command.full <- paste(bs.command_start, 'list project -f csv | tr "," "\t" 1>>bs_data.tsv')
  runExternal(command.full)
  outputnames[1, ] <- c("bs_data.tsv", "basespace-projects.tsv")
}

if (action == "dir") {
  command.full <- paste(bs.command_start, "dir dataset -f csv --name", name, '  | tr "," "\t" 1>>bs_data.tsv')
  runExternal(command.full)
  outputnames[1, ] <- c("bs_data.tsv", paste(name, "-dir.tsv", sep=""))
}

if (action == "info") {
  command.full <- paste(bs.command_start, "get dataset --name", name, " 1>>bs_data.txt")
  runExternal(command.full)
  outputnames[1, ] <- c("bs_data.txt", paste(name, "-info.txt", sep=""))
}

if (action == "download_dataset") {
  command.full <- paste(bs.command_start, "download dataset --name", name, " -z -o bs_download")
  runExternal(command.full)
  outputnames[1, ] <- c("bs_data.tar", paste(name, ".tar", sep=""))
}

if (action == "download_project") {
  command.full <- paste(bs.command_start, "download project --name", name, " --extension=fastq.gz -z -o bs_download")
  runExternal(command.full)
  outputnames[1, ] <- c("bs_data.tar", paste(name, ".tar", sep=""))
}

# if (action == "download_dataset_id") {
#   command.full <- paste(bs.command_start, "download dataset -i", name, " -z -o bs_download")
#   runExternal(command.full)
# }

system("ls -lah")


if (file.exists("bs_download.tar.gz")) {

  runExternal("mkdir extracted_files")
  # extract without path
  print("extract bs_download.tar.gz")
  runExternal("tar xzf bs_download.tar.gz --transform='s/.*\\///' -C extracted_files")

  system("ls -lah extracted_files")

  # add only .fastq.gz files (no .json) to a tar package (without compression, because individual files are compressed already)
  # use bash to fill in the wildcard '*'
  print("package fastq files")

  fastq.files <- list.files(path = "extracted_files", pattern=".fastq.gz$")
  if (length(fastq.files) == 0) {
    stop("No FASTQ files in the project.")
  } else {
    joined_files <- paste(fastq.files, collapse=" ")
    cmd <- paste("tar cf bs_data.tar -C extracted_files", joined_files)
    runExternal(cmd)
  }
}

# Write output definitions file
write_output_definitions(outputnames)
