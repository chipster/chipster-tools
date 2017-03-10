# TOOL rseqc_infer_rnaseq_experiment.R: "RNA-seq strandedness inference and inner distance estimation using RseQC" (Given FASTQ files, this tool aligns a subset of the reads against a reference genome. Alignments are then compared to reference annotation to infer strandedness. Please see the manual for help with interpreting the results. For paired-end reads the inner distance distribution is also calculated.)
# INPUT reads1.fq: "Read 1 FASTQ" TYPE GENERIC
# INPUT OPTIONAL reads2.fq: "Read 2 FASTQ" TYPE GENERIC
# OUTPUT experiment_data.txt
# OUTPUT OPTIONAL inner_distance.pdf
# PARAMETER organism: "Genome" TYPE ["FILES genomes/bed .bed"] DEFAULT "SYMLINK_TARGET genomes/bed/default .bed" (Genome or transcriptome that you would like to align your reads against.)

# check out if the file is compressed and if so unzip it
source(file.path(chipster.common.path, "zip-utils.R"))
unzipIfGZipFile("reads1.fq")
unzipIfGZipFile("reads2.fq")

# Is the submitted data paired end?
pe <- FALSE
if (file.exists("reads2.fq")){
	pe <- TRUE
}

# Make subsets of FASTQ files for faster processing
seqtk.binary <- file.path(chipster.tools.path, "seqtk", "seqtk")
system(paste(seqtk.binary, "sample -s 15 reads1.fq 200000 > subset.reads1.fq"))
if (pe){
    system(paste(seqtk.binary, "sample -s 15 reads2.fq 200000 > subset.reads2.fq"))
}

# Align against reference genome
bowtie.binary <- c(file.path(chipster.tools.path, "bowtie2", "bowtie2"))
bowtie.genome <- c(file.path(chipster.tools.path, "genomes", "indexes", "bowtie2", organism))
bowtie.command <- paste("bash -c '", bowtie.binary, "-p", chipster.threads.max, "-x", bowtie.genome)
if (pe){
	bowtie.command <- paste(bowtie.command, "-1 subset.reads1.fq -2 subset.reads2.fq")
}else{
	bowtie.command <- paste(bowtie.command, "-U subset.reads1.fq")
}
bowtie.command <- paste(bowtie.command, "-S alignment.sam 2>> bowtie2.log'")
system(bowtie.command)

internal.bed <- file.path(chipster.tools.path, "genomes", "bed", paste(organism, ".bed" ,sep="" ,collapse=""))

# Infer experiment
# tools-bin RSeQC
ie.binary <- c(file.path(chipster.tools.path, "rseqc", "infer_experiment.py"))

ie.command <- paste(ie.binary, "-i alignment.sam -r", internal.bed, "> experiment_data.txt")
system(ie.command)

# Add some helpful explanation to the output
# Get the values from the output file
system(paste(" grep \":\" experiment_data.txt | awk -F \": \" '{print $2}' > values.txt"))
# read them into R
values <- scan("values.txt", what=double(), sep="\n")
# Make sure we dont get a division by zero error
values[2] <- values[2] + 0.0001
values[3] <- values[3] + 0.0001
# Check the case
ratio <- values[2] / values[3]


# Add additional help lines to output
if (pe){
	if (ratio > 10){
		message <- paste("\nIt seems the data is stranded. Read 1 is always on the same strand as the gene.")
		message <- paste(message, "\n\nCorresponding parameters are:")
		message <- paste(message, "\nTopHat, Cufflinks and Cuffdiff: library-type fr-secondstrand")
		message <- paste(message, "\nHTSeq: stranded -- yes")
		message <- paste(message, "\nRSeQC: 1++,1–,2+-,2-+" )
	}else if (ratio < 0.1){
		message <- paste("\nIt seems the data is stranded. Read 2 is always on the same strand as the gene.")
		message <- paste(message, "\n\nCorresponding parameters are:")
		message <- paste(message, "\nTopHat, Cufflinks and Cuffdiff: library-type fr-firststrand")
		message <- paste(message, "\nHTSeq: stranded -- reverse")
		message <- paste(message, "\nRSeQC: 1+-,1-+,2++,2–" )
	}else{
		message <- paste("\nIt seems the data is unstranded.")
		message <- paste(message, "\n\nCorresponding parameters are:")
		message <- paste(message, "\nTopHat, Cufflinks and Cuffdiff: library-type fr-unstranded")
		message <- paste(message, "\nHTSeq: stranded -- no")
		message <- paste(message, "\nRSeQC: none")
	}	
}else{
	if (ratio > 10){
		message <- paste("\nIt seems the data is stranded. Read is always on the same strand as the gene.")
		message <- paste(message, "\n\nCorresponding parameters are:")
		message <- paste(message, "\nTopHat, Cufflinks and Cuffdiff: library-type fr-secondstrand")
		message <- paste(message, "\nHTSeq: stranded -- yes")
		message <- paste(message, "\nRSeQC: ++,–" )
	}else if (ratio < 0.1){
		message <- paste("\nIt seems the data is stranded. Read is always on the opposite strand to gene.")
		message <- paste(message, "\n\nCorresponding parameters are:")
		message <- paste(message, "\nTopHat, Cufflinks and Cuffdiff: library-type fr-firststrand")
		message <- paste(message, "\nHTSeq: stranded -- reverse")
		message <- paste(message, "\nRSeQC: +-,-+" )
	}else{
		message <- paste("\nIt seems the data is unstranded.")
		message <- paste(message, "\n\nCorresponding parameters are:")
		message <- paste(message, "\nTopHat, Cufflinks and Cuffdiff: library-type fr-unstranded")
		message <- paste(message, "\nHTSeq: stranded -- no")
		message <- paste(message, "\nRSeQC: none")
	}	
}	
# Add input names to output
source(file.path(chipster.common.path, "tool-utils.R"))
inputnames <- read_input_definitions()

message <- paste(message, "\n\nInput files were assigned as follows:")
message <- paste(message, "\nRead 1 file:", inputnames$reads1.fq)
if(pe){
	message <- paste(message, "\nRead 2 file:", inputnames$reads2.fq)
}

write(message, file = "experiment_data.txt", append = TRUE)

# Inner distance. Only for paired end reads
if (pe){
    # tools-bin RSeQC
    ie.binary <- c(file.path(chipster.tools.path, "rseqc", "inner_distance.py"))
   
    id.command <- paste(ie.binary, "-i alignment.sam -r", internal.bed, " -o id")
    system(id.command)
    try(source("id.inner_distance_plot.r"), silent=TRUE)
    system("mv id.inner_distance_plot.pdf inner_distance.pdf")
}

