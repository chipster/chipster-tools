# TOOL bwa-mem-paired-end.R: "BWA MEM for paired-end reads" (Aligns reads to genomes using the BWA MEM algorithm. Results are sorted and indexed BAM files.
# Note that this BWA tool uses publicly available genomes. If you would like to align reads against your own datasets, please use the tool \"BWA MEM for single end reads and own genome\".)
# INPUT reads1.txt: "Paired-end read set 1 to align" TYPE GENERIC
# INPUT reads2.txt: "Paired-end read set 2 to align" TYPE GENERIC
# OUTPUT bwa.bam
# OUTPUT bwa.log
# OUTPUT OPTIONAL bwa.bam.bai
# PARAMETER organism: "Genome or transcriptome" TYPE [Arabidopsis_thaliana.TAIR10, Bos_taurus.UMD3.1, Canis_familiaris.CanFam3.1, Drosophila_melanogaster.BDGP6, Felis_catus.Felis_catus_6.2, Gallus_gallus.Galgal4, Gallus_gallus.Gallus_gallus-5.0, Gasterosteus_aculeatus.BROADS1, Halorubrum_lacusprofundi_atcc_49239.ASM2220v1, Homo_sapiens.GRCh37.75, Homo_sapiens.GRCh38, Homo_sapiens_mirna, Medicago_truncatula.MedtrA17_4.0, Mus_musculus.GRCm38, Mus_musculus_mirna, Oryza_sativa.IRGSP-1.0, Ovis_aries.Oar_v3.1, Populus_trichocarpa.JGI2.0, Rattus_norvegicus_mirna, Rattus_norvegicus.Rnor_5.0, Rattus_norvegicus.Rnor_6.0, Schizosaccharomyces_pombe.ASM294v2, Solanum_tuberosum.SolTub_3.0, Sus_scrofa.Sscrofa10.2, Vitis_vinifera.IGGP_12x, Yersinia_enterocolitica_subsp_palearctica_y11.ASM25317v1, Yersinia_pseudotuberculosis_ip_32953_gca_000834295.ASM83429v1] DEFAULT Homo_sapiens.GRCh38 (Genome or transcriptome that you would like to align your reads against.)
# PARAMETER OPTIONAL index.file: "Create index file" TYPE [index_file: "Create index file", no_index: "No index file"] DEFAULT no_index (Creates index file for BAM. By default no index file.)
# PARAMETER mode: "Data source" TYPE [normal: " Illumina, 454, IonTorrent reads longer than 70 base pairs", pacbio: "PacBio subreads"] DEFAULT normal (Defining the type of reads will instruct the tool to use a predefined set of parameters optimized for that read type.)
# RUNTIME R-4.1.1

# KM 11.11.2014

# check out if the file is compressed and if so unzip it
source(file.path(chipster.common.lib.path, "zip-utils.R"))
unzipIfGZipFile("reads1.txt")
unzipIfGZipFile("reads2.txt")

# bwa
bwa.binary <- file.path(chipster.tools.path, "bwa", "bwa mem")
bwa.genome <- file.path(chipster.tools.path, "genomes", "indexes", "bwa", organism)
command.start <- paste("bash -c '", bwa.binary)

mode.parameters <- ifelse(mode == "pacbio", "-x pacbio", "")

# command ending
command.end <- paste(bwa.genome, "reads1.txt reads2.txt 1> alignment.sam 2>> bwa.log'")

# run bwa alignment
bwa.command <- paste(command.start, mode.parameters, command.end)

echo.command <- paste("echo '", bwa.binary, mode.parameters, bwa.genome, "reads1.txt reads2.txt' > bwa.log")
# stop(paste('CHIPSTER-NOTE: ', bwa.command))
runExternal(echo.command)
runExternal(bwa.command)

# samtools binary
samtools.binary <- c(file.path(chipster.tools.path, "samtools", "bin", "samtools"))

# convert sam to bam
runExternal(paste(samtools.binary, "view -bS alignment.sam -o alignment.bam"))

# sort bam
runExternal(paste(samtools.binary, "sort alignment.bam -o alignment.sorted.bam"))

# index bam
syrunExternalstem(paste(samtools.binary, "index alignment.sorted.bam"))

# rename result files
runExternal("mv alignment.sorted.bam bwa.bam")
if (index.file == "index_file") {
  runExternal("mv alignment.sorted.bam.bai bwa.bam.bai")
}

# Handle output names
#
source(file.path(chipster.common.lib.path, "tool-utils.R"))

# read input names
inputnames <- read_input_definitions()

# Determine base name
base1 <- strip_name(inputnames$reads1.txt)
base2 <- strip_name(inputnames$reads2.txt)
basename <- paired_name(base1, base2)

# Make a matrix of output names
outputnames <- matrix(NA, nrow = 2, ncol = 2)
outputnames[1, ] <- c("bwa.bam", paste(basename, ".bam", sep = ""))
outputnames[2, ] <- c("bwa.bam.bai", paste(basename, ".bam.bai", sep = ""))

# Write output definitions file
write_output_definitions(outputnames)
