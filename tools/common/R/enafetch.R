# TOOL enafetch.R: "Retrieve datasets from ENA database" (Retrieve data from the ENA database based on the entry ID or name.)
# OUTPUT OPTIONAL enafetch.log
# OUTPUT OPTIONAL {...}.gz: "gzipped files"
# OUTPUT OPTIONAL {...}.fasta: "fasta files"
# OUTPUT OPTIONAL {...}.dat: "embl files"
# OUTPUT OPTIONAL {...}.txt: "text files"
# OUTPUT OPTIONAL {...}.xml: "xml files"
# OUTPUT OPTIONAL {...}.bam: "bam files"
# OUTPUT OPTIONAL {...}.cram: "cram files"
# OUTPUT OPTIONAL {...}.crai: "cram files"
# PARAMETER entry_id: "ENA ID" TYPE STRING DEFAULT "entry" (Give the ID of the ENA dataset to be retrieved. For example: SRR000021) 
# PARAMETER format: "Data to retrieve" TYPE [default: "Default format", fastq: "fastq reads", fasta: "fasta sequence", embl: "embl formatted sequence" ] DEFAULT default (Define the data forma for retried data. Not that if you retrieve NGS reads, you should use fastq format in stead of default format)
# PARAMETER index: "Download CRAM index" TYPE [yes: "Yes", no: "No" ] DEFAULT no ( Download CRAM index files with submitted CRAM files, if any. This selection is ignored if fastq format is selected.)
# PARAMETER OPTIONAL save_log: "Collect a log file" TYPE [yes: yes, no: no] DEFAULT no (Collect a log file.)

# KM 21.03.2019


ena.path <- file.path(chipster.tools.path, "enabrowsertools/python3")
#ena.path <- ("/opt/chipster/tools_local/enaBrowserTools-1.5.4/python")
#turn of cacheing
ena.binary <- file.path(ena.path, "enaDataGet")

index.option <- (" ")

if ( index == "yes"){
	if ( format != "fastq" ){
    	index.option <- ("-i")
	}
}

if ( format != "default" ){
command.full <- paste(ena.binary, index.option, '-f', format, entry_id, ' 1>>enafetch.log 2>>enafetch.log')
}else{
	command.full <- paste(ena.binary, index.option, entry_id, ' 1>>enafetch.log 2>>enafetch.log')
}
cat(command.full, "\n", file="enafetch.log", append=TRUE)
system(command.full)
#temporary bug fix. Command needs to be run twice to get both fastq files
#system(command.full)

# fasta and dat files are downloaded to the current directory
# gunzip them if needed 
system("gunzip *.gz")


##check storage resources and needs
#freespace <- system(" df ./ | awk '{print $4}' ", intern = TRUE )
#sra.binary <- file.path(sra.path, "vdb-dumb")
#commad.full <- paste(sra.binary, ' --info ', entry_id, '| grep -w size | tr -d ","', " | awk '{print $3}'")
#srasize <-system(command.full)
#space_needed = 4*srasize

#if ( freespace < space_needed ){
#	stop("CHIPSTER-NOTE: THe SRA entry you are trying to retrieve is too larage for the chipster server")
#}

#fastq files are stored to a subdrectory. Copy them to current directory (but keep them gzipped)
system("cp ./*/* ./")
system("cp ./*/*/* ./")
## check that we got something
#gz.files <- Sys.glob("*.gz")
#if (length(gz.files) == 0) {
#	system('echo enafetch.log:')
#	system('cat enafetch.log')
#	stop("Fetching fastq files failed, see the end of the details for more information")
#}

system("ls -l >> enafetch.log")

## rename outputs and create output names file
#source(file.path(chipster.common.path, "tool-utils.R"))
#outputnames <- matrix(NA, nrow=length(gz.files), ncol=2)
#for (i in 1:length(gz.files)) {
#	original.name <- gz.files[i]
#	new.name <- paste("ena_", i, ".gz", sep ="")
#	outputnames[i,] <- c(new.name, original.name)
#	system(paste("mv", original.name, new.name))
#}
#write_output_definitions(outputnames)

cat("--------------------", "\n", file="enafetch.log", append=TRUE)
system("ls -l >> enafetch.log")

if ( save_log == "no") {
	system ("rm -f enafetch.log")
}
