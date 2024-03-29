# TOOL prinseq-statistics.R: "Read quality statistics with PRINSEQ" (Calculates general statistics of the reads in the given FASTQ file. This tool is based on the PRINSEQ program. Please note that if your file is larger than 4 GB, we recommend that you submit only a sample of reads for the quality statistics analysis, because PRINSEQ uses a lot of memory when producing the html report and might fail with bigger files. You can use the tool Utilities / Make a subset of FASTQ for this.)
# INPUT fastqfile: "Input reads file" TYPE GENERIC
# OUTPUT OPTIONAL reads-stats.html
# PARAMETER OPTIONAL input.mode: "Input file format" TYPE [ fq: "FASTQ", fa: "FASTA"] DEFAULT fq (Define the file format of the reads file)

# KM 17.1.2012
# EK 23.4.2012
# AMS 17.2.2014, removed table format output, added graph_stats parameter to cl

# check out if the file is compressed and if so unzip it
source(file.path(chipster.common.lib.path, "zip-utils.R"))
unzipIfGZipFile("fastqfile")

# binary
binary.stats <- c(file.path(chipster.tools.path, "prinseq", "prinseq-lite.pl"))

# quality encoding check
# quality.scale <- ifelse(phred64 == "y", "-phred64", "")

# command to generate result table
# system('printf "%s\t%s\t%s\n"  Class Feature Value > reads-stats.tsv')
# if (input.mode == "fq") {
# 	command.stats <- paste("perl", binary.stats, " -fastq fastqfile -out_good null -out_bad null -stats_all >> reads-stats.tsv")
# }
# if (input.mode == "fa") {
# 	command.stats <- paste("perl", binary.stats, " -fasta fastqfile -out_good null -out_bad null -stats_all >> reads-stats.tsv")
# }
#
# ret <- system(command.stats)
# if (ret > 0) {
# 	stop('Unsupported input file type, please see tool output for more details.')
# }

# commands to generate graph file
if (input.mode == "fq") {
   command.graph <- paste("perl", binary.stats, " -fastq fastqfile -out_good null -out_bad null -graph_stats ld,gc,qd,ns,pt,ts,de -graph_data tmp_graph_file ")
}
if (input.mode == "fa") {
   command.graph <- paste("perl", binary.stats, " -fasta fastqfile -out_good null -out_bad null -graph_stats ld,gc,qd,ns,pt,ts,de -graph_data tmp_graph_file ")
}
system(command.graph)

# create html file
binary.graph <- c(file.path(chipster.tools.path, "prinseq", "prinseq-graphs.pl"))
command.graph <- paste("perl", binary.graph, " -i tmp_graph_file -html_all -o reads-stats")
system(command.graph)

# Handle output names
source(file.path(chipster.common.lib.path, "tool-utils.R"))

# read input names
inputnames <- read_input_definitions()

basename <- strip_name(inputnames$fastqfile)

# Make a matrix of output names
outputnames <- matrix(NA, nrow = 1, ncol = 2)
outputnames[1, ] <- c("reads-stats.html", paste(basename, "_prinseq.html", sep = ""))

# Write output definitions file
write_output_definitions(outputnames)
