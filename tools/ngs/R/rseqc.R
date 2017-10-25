# TOOL rseqc.R: "RNA-seq quality metrics with RseQC" (Given an RNA-seq BAM file and gene and exon locations in a BED file, this tool reports several quality metrics such as coverage uniformity, gene and junction saturation, junction annotation and alignment statistics. You can provide your own BED file or use one of the internal annotations. This tool is based on the RSeQC package.)
# INPUT alignment.bam: "BAM file" TYPE GENERIC
# INPUT OPTIONAL reference_file: "BED file" TYPE GENERIC
# OUTPUT OPTIONAL RSeQC.txt
# OUTPUT OPTIONAL RSeQC_report.pdf
# PARAMETER organism: "Organism" TYPE [other: "Own BED file", "FILES genomes/bed .bed"] DEFAULT other (Choose one of the reference organisms or provide your own BED file.)
# PARAMETER chr: "Chromosome names in my BAM file look like" TYPE [chr1, 1] DEFAULT 1 (Chromosome names must match in the BAM file and in the annotation. Check your BAM and choose accordingly. This only applies if you have generated your BAM outside Chipster and are using one of the reference organisms.)
# PARAMETER OPTIONAL rpkm: "Generate RPKM saturation plot" TYPE [yes, no] DEFAULT no (BAM file containing more than 100 million alignments will make this analysis very slow. Try disabling it if RSeQC takes a very long time or fails to complete.)
# PARAMETER OPTIONAL paired: "Generate inner distance plot" TYPE [yes, no] DEFAULT no (Calculate the inner distance (or insert size\) between two paired RNA reads. The distance is the mRNA length between two paired fragments.)


# PARAMETER OPTIONAL strandedness: "Strandedness of reads" TYPE [none, 1++1--,2+-2-+: "1++,1--,2+-,2-+", 1+-1-+2++2--: "1+-,1-+,2++,2--", ++--: "++,--", +--+: "+-,-+"] DEFAULT none (How reads were stranded during sequencing. This only needed for RPKM saturation plot. If unsure, use tool RNA-seq strandedness inference and inner distance estimation using RseQC.)


# AMS 09.01.2014
# AMS 23.05.2014 added inner distance plot
# AMS 03.12.2014 improved error handling for the plots
# AMS 07.04.2015, combined pdf outputs
# AMS 21.11.2016, tools-bin Python, RSeQC 2.6.4, internal BED files

rseqc.path <- c(file.path(chipster.tools.path, "rseqc"))

samtools.binary <- c(file.path(chipster.tools.path, "samtools", "samtools"))
system(paste(samtools.binary, "index alignment.bam > alignment.bam.bai"))

# Which BED to use.
# If user has provided a BED, we use it.
if (file.exists("reference_file")){
	bed.file <- paste("reference_file")
}else{
	if (organism == "other"){
		stop(paste('CHIPSTER-NOTE: ', "Choose one of the reference organisms or provide your own BED file."))
	}else{
		# If not, we use the internal one.
		internal.bed <- file.path(chipster.tools.path, "genomes", "bed", paste(organism, ".bed" ,sep="" ,collapse=""))
		# If chromosome names in BAM have chr, we make a temporary copy of BED with chr names, otherwise we use it as is.
		if(chr == "chr1"){
			source(file.path(chipster.common.path, "gtf-utils.R"))
			addChrToGtf(internal.bed, "internal_chr.bed") 
			bed.file <- paste("internal_chr.bed")
		}else{
			bed.file <- paste(internal.bed)
		}		
	}
}


# geneBody_coverage
binary <- c(file.path(rseqc.path, "geneBody_coverage.py"))
command <- paste(binary, "-i alignment.bam -r", bed.file, "-o RSeQC")
system(command)
try(source("RSeQC.geneBodyCoverage.r"), silent=TRUE)
# Outputs are renamed to control their order in the joined PDF
system("mv RSeQC.geneBodyCoverage.curves.pdf 01.pdf")

# junction_saturation
binary <- c(file.path(rseqc.path, "junction_saturation.py"))
command <- paste(binary, "-i alignment.bam -r", bed.file, "-o RSeQC")
system(command)
try(source("RSeQC.junctionSaturation_plot.r"), silent=TRUE)
# Outputs are renamed to control their order in the joined PDF
system("mv RSeQC.junctionSaturation_plot.pdf 02.pdf")

# junction_annotation
binary <- c(file.path(rseqc.path, "junction_annotation.py"))
command <- paste(binary, "-i alignment.bam -r", bed.file, "-o RSeQC")
system(command)
try(source("RSeQC.junction_plot.r"), silent=TRUE)
# Outputs are renamed to control their order in the joined PDF
system("mv RSeQC.splice_events.pdf 03.pdf")
system("mv RSeQC.splice_junction.pdf 04.pdf")

#RPKM_saturation
binary <- c(file.path(rseqc.path, "RPKM_saturation.py"))
command <- paste(binary, "-i alignment.bam -r", bed.file, "-o RSeQC")
system(command)
try(source("RSeQC.saturation.r"), silent=TRUE)
# Outputs are renamed to control their order in the joined PDF
system("mv RSeQC.saturation.pdf 05.pdf")

# bam_stat
system("echo 'bam_stats:' > RSeQC.txt")
system("echo >> RSeQC.txt")
binary <- c(file.path(rseqc.path, "bam_stat.py"))
command <- paste(binary, "-i alignment.bam >> RSeQC.txt")
system(command)

# read_distribution
system("echo \"\n\n\" >> RSeQC.txt")
system("echo 'read_distribution:' >> RSeQC.txt")
system("echo >> RSeQC.txt")
binary <- c(file.path(rseqc.path, "read_distribution.py"))
command <- paste(binary, "-i alignment.bam -r", bed.file, ">> RSeQC.txt")
system(command)

# inner_distance
if (paired == "yes"){
	binary <- c(file.path(rseqc.path, "inner_distance.py"))
	command <- paste(binary, "-i alignment.bam -r", bed.file, "-o RSeQC")
	system(command)
	try(source("RSeQC.inner_distance_plot.r"), silent=TRUE)
	# Outputs are renamed to control their order in the joined PDF
	system("mv RSeQC.inner_distance_plot.pdf 06.pdf")
}

# Join the PDFs 
system("gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=RSeQC_report.pdf *.pdf")

# Handle output names
source(file.path(chipster.common.path, "tool-utils.R"))

# read input names
inputnames <- read_input_definitions()

base <- strip_name(inputnames$alignment.bam)

# Make a matrix of output names
outputnames <- matrix(NA, nrow=2, ncol=2)
outputnames[1,] <- c("RSeQC.txt", paste(base, "_rseqc.txt", sep=""))
outputnames[2,] <- c("RSeQC_report.pdf", paste(base, "_rseqc.pdf", sep=""))

# Write output definitions file
write_output_definitions(outputnames)
