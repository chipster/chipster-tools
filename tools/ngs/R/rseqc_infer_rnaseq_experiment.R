# TOOL rseqc_infer_rnaseq_experiment.R: "Infer strandedness and inner distance from FASTQ" (Subsets of the input FASTQ files are first aligned against reference genomes. Alignments are then compared to reference annotation to infer strandedness. Please see the manual for help with interpreting the results. For paired-end reads matepair inner distance distribution is also calculated.)
# INPUT reads1.fq: "Read 1 FASTQ" TYPE GENERIC
# INPUT OPTIONAL reads2.fq: "Read 2 FASTQ" TYPE GENERIC
# OUTPUT experiment_data.txt
# OUTPUT OPTIONAL inner_distance.pdf
# PARAMETER organism: "Genome" TYPE [Arabidopsis_thaliana.TAIR10.32, Bos_taurus.UMD3.1.86, Canis_familiaris.BROADD2.67, Canis_familiaris.CanFam3.1.86, Drosophila_melanogaster.BDGP5.78, Drosophila_melanogaster.BDGP6.86, Felis_catus.Felis_catus_6.2.86, Gallus_gallus.Galgal4.85, Gallus_gallus.Gallus_gallus-5.0.86, Gasterosteus_aculeatus.BROADS1.86, Halorubrum_lacusprofundi_atcc_49239.ASM2220v1.32, Homo_sapiens.GRCh37.75, Homo_sapiens.GRCh38.86, Homo_sapiens.NCBI36.54, Medicago_truncatula.MedtrA17_4.0.32, Mus_musculus.GRCm38.86, Mus_musculus.NCBIM37.67, Oryza_sativa.IRGSP-1.0.32, Ovis_aries.Oar_v3.1.86, Populus_trichocarpa.JGI2.0.32, Rattus_norvegicus.RGSC3.4.69, Rattus_norvegicus.Rnor_5.0.79, Rattus_norvegicus.Rnor_6.0.86, Schizosaccharomyces_pombe.ASM294v2.32, Solanum_tuberosum.SolTub_3.0.32, Sus_scrofa.Sscrofa10.2.86, Vitis_vinifera.IGGP_12x.32, Yersinia_enterocolitica_subsp_palearctica_y11.ASM25317v1.32, Yersinia_pseudotuberculosis_ip_32953_gca_000834295.ASM83429v1.32] DEFAULT Homo_sapiens.GRCh38.86 (Genome or transcriptome that you would like to align your reads against.)

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
#stop("CHIPSTER-NOTE:", bowtie.command)
system(bowtie.command)

internal.bed <- file.path(chipster.tools.path, "genomes", "bed", paste(organism, ".bed" ,sep="" ,collapse=""))

# Infer experiment
# tools-bin RSeQC
ie.binary <- c(file.path(chipster.tools.path, "rseqc", "infer_experiment.py"))
# Old RSeQC
#ie.binary <- c(file.path(chipster.tools.path, "RSeQC", "scripts", "infer_experiment.py"))

ie.command <- paste(ie.binary, "-i alignment.sam -r", internal.bed, "> experiment_data.txt")
system(ie.command)

# Add some helpful explanation to the output
# Get the values from the output file
system(paste(" grep 1++ experiment_data.txt | awk -F \": \" '{print $2}' > values.txt"))
system(paste(" grep 2++ experiment_data.txt | awk -F \": \" '{print $2}' >> values.txt"))
# read them into R
values <- scan("values.txt", what=double(), sep="\n")
# Make sure we dont get a division by zero error
values[1] <- values[1] + 0.0001
values[2] <- values[2] + 0.0001
# Check the case
ratio <- values[1] / values[2]

if (ratio > 10){
	message <- paste("\nIt seems the data is stranded. Read 1 is always on the same strand as the gene.")
}else if (ratio < 0.1){
	message <- paste("\nIt seems the data is stranded. Read 2 is always on the same strand as the gene.")
}else{
	message <- paste("\nIt seems the data is unstranded.")	
}
write(message, file = "experiment_data.txt", append = TRUE)

# Inner distance. Only for paired end reads
if (pe){
    # tools-bin RSeQC
    ie.binary <- c(file.path(chipster.tools.path, "rseqc", "inner_distance.py"))
    # Old RSeQC
    #ie.binary <- c(file.path(chipster.tools.path, "RSeQC", "scripts", "inner_distance.py"))

    id.command <- paste(ie.binary, "-i alignment.sam -r", internal.bed, " -o id")
    system(id.command)
    try(source("id.inner_distance_plot.r"), silent=TRUE)
    system("mv id.inner_distance_plot.pdf inner_distance.pdf")
}

