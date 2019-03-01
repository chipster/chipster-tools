# TOOL illumina_download.R: "Retrieve data from Illumina BaseSpace" (Retrieve data from illumina BaseSpace. This tool requires that you have an access token for bs clinet)
# OUTPUT OPTIONAL bs.log
# OUTPUT OPTIONAL bs_data.tsv
# OUTPUT OPTIONAL bs_data.txt
# OUTPUT OPTIONAL {...}.fastqc.gz: "FASTQ files"
# OUTPUT OPTIONAL bs_download.tar.gz
# PARAMETER OPTIONAL name: "Name or dataset" TYPE STRING (Give the name of the dataset. This paramter is not needed if you use this tool to just list your datasets in Illumina BaseSpace) 
# PARAMETER action: "Action" TYPE [list: "List datasets", download: "Download dataset", dir: "Display content of a dataset", info: "Display detailed information about a dataset"  ] DEFAULT list (Action to be performed)
# PARAMETER apiserver: "Api-server" TYPE [api.basespace.illumina.com: "api.basespace.illumina.com"] DEFAULT api.basespace.illumina.com (Define the BaseSpace server to be used) 
# PARAMETER token: "Access token" TYPE STRING (Your personal Illumina BaceSpace access token) 
 
# KM 25.02.2018


bs.binary <- file.path("/opt/chipster/tools_local", "bs/bin/bs")
#turn of cacheing
bs.command_start <- paste(bs.binary, " --api-server=https://", apiserver, " --access-token=", token, sep="" )

if (action== "list"){
  command.full <- paste(bs.command_start, 'list dataset -f csv | tr "," "\t" 1>>bs_data.tsv 2>>bs.log')
  cat(command.full, "\n", file="bs.log", append=TRUE)
  system(command.full)
}

if (action== "dir"){
	command.full <- paste(bs.command_start, 'dir dataset -f csv --name', name, '  | tr "," "\t" 1>>bs_data.tsv 2>>bs.log')
	cat(command.full, "\n", file="bs.log", append=TRUE)
	system(command.full)
}

if (action== "info"){
	command.full <- paste(bs.command_start, 'get dataset --name', name, '  1>>bs_data.txt 2>>bs.log')
	cat(command.full, "\n", file="bs.log", append=TRUE)
	system(command.full)
}

if (action== "download"){
	command.full <- paste(bs.command_start, 'download dataset -o download --name', name, ' -z -o bs_download 1>>bs.log 2>>bs.log')
	cat(command.full, "\n", file="bs.log", append=TRUE)
	system(command.full)
}

system("ls -l >> bs.log")

# check that we got something
#"fastq.files <- Sys.glob("*.fastq.gz")
#if (length(fastq.files) == 0) {
#	system('echo srafetch.log:')
#	system('cat srafetch.log')
#	stop("Fetching fastq files failed, see the end of the details for more information")
#}

## rename outputs and create output names file
#source(file.path(chipster.common.path, "tool-utils.R"))
#outputnames <- matrix(NA, nrow=length(fastq.files), ncol=2)
#for (i in 1:length(fastq.files)) {
#	original.name <- fastq.files[i]
#	new.name <- paste("sra_reads_", i, ".fastqc.gz", sep ="")
#	outputnames[i,] <- c(new.name, original.name)#
#	system(paste("mv", original.name, new.name))
#}
#write_output_definitions(outputnames)

