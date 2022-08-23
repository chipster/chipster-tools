# TOOL star-paired-end.R: "STAR for paired end reads" (Aligns paired end RNA-seq reads to a genome. If you have just one pair of read files, Chipster sets reads 1 file and reads 2 file based on file names. If you have more pairs of read files for one sample, you need to provide a list of filenames of the FASTQ files for each direction \(e.g. 1files.txt and 2files.txt\). You can generate the lists with the tool \"Utilities \\\ Make a list of filenames\". Alignment results are given in a BAM file, which is automatically indexed and hence ready to be viewed in Chipster genome browser.)
# INPUT reads{...}.fq: "Reads" TYPE GENERIC
# INPUT OPTIONAL reads1.txt: "List of read 1 files" TYPE GENERIC
# INPUT OPTIONAL reads2.txt: "List of read 2 files" TYPE GENERIC
# INPUT OPTIONAL annotation.gtf: "Optional GTF file" TYPE GENERIC
# OUTPUT OPTIONAL alignment.bam
# OUTPUT OPTIONAL alignment.bam.bai
# OUTPUT OPTIONAL Log_progress.txt
# OUTPUT OPTIONAL Log_final.txt
# PARAMETER organism: "Genome" TYPE [Homo_sapiens.GRCh38.95, Mus_musculus.GRCm38.95, Rattus_norvegicus.Rnor_6.0.95] DEFAULT Homo_sapiens.GRCh38.95 (Genome that you would like to align your reads against.)
# PARAMETER OPTIONAL index.file: "Create index file" TYPE [index_file: "Create index file", no_index: "No index file"] DEFAULT no_index (Creates index file for BAM. By default no index file.)
# PARAMETER OPTIONAL alignments.per.read: "Maximum alignments per read" TYPE INTEGER DEFAULT 10 (Maximum number of multiple alignments allowed for a read: if exceeded, the read is considered unmapped.)
# PARAMETER OPTIONAL mismatches.per.pair: "Maximum mismatches per alignment" TYPE INTEGER DEFAULT 10 (Maximum number of mismatches per alignment. Use value 999 to switch off this filter.)
# PARAMETER OPTIONAL out.filter.mismatch.nover.lmax: "Mismatch ratio" TYPE DECIMAL DEFAULT 0.3 (Alignment will be output only if its ratio of mismatches to mapped length is less than or equal to this value.)
# PARAMETER OPTIONAL align.intron.min: "Minimum intron size" TYPE INTEGER DEFAULT 21 (Minimum intron size.)
# PARAMETER OPTIONAL align.intron.max: "Maximum intron size" TYPE INTEGER DEFAULT 0 (If 0, max intron size will be determined automatically, please see the manual page.)
# PARAMETER OPTIONAL align.mates.gap.max: "Maximum gap between two mates" TYPE INTEGER DEFAULT 0 (If 0, max intron gap will be determined automatically, please see the manual page.)
# PARAMETER OPTIONAL log.files: "Create log files" TYPE [final_log: "Final log only", final_and_progress: "Final and progress logs", no_logs: "No logs"] DEFAULT final_log (Do you want to create a log file? By default only the final log is created.)
# SLOTS 5

source(file.path(chipster.common.path,"tool-utils.R"))
source(file.path(chipster.common.path,"bam-utils.R"))
source(file.path(chipster.common.path,"zip-utils.R"))

# check out if the file is compressed and if so unzip it
input.names <- read.table("chipster-inputs.tsv",header = FALSE,sep = "\t")
for (i in 1:nrow(input.names)) {
  unzipIfGZipFile(input.names[i,1])
}

# setting up STAR
# latest STAR is not compatible with old indexes
# use the latest version, set runtime R-4.1.1 and latest samtools after indexes are updated
# star.binary <- c(file.path(chipster.tools.path,"STAR","STAR"))
star.binary <- c(file.path(chipster.tools.path,"STAR-2.5.3a","STAR"))
path.star.index <- c(file.path(chipster.tools.path,"genomes","indexes","star",organism))
path.gtf <- c(file.path(chipster.tools.path,"genomes","gtf",organism))
samtools.binary <- c(file.path(chipster.tools.path,"samtools-0.1.19","samtools"))
#samtools.binary <- c(file.path(chipster.tools.path, "samtools", "bin", "samtools"))

# Input files
if (fileOk("reads1.txt",0) && fileOk("reads2.txt",0)) {
  # Case: list files exist
  reads1.list <- make_input_list("reads1.txt")
  reads2.list <- make_input_list("reads2.txt")
  if (identical(intersect(reads1.list,reads2.list),character(0))) {
    reads1 <- paste(reads1.list,sep = "",collapse = ",")
    reads2 <- paste(reads2.list,sep = "",collapse = ",")
  } else {
    stop(paste('CHIPSTER-NOTE: ',"One or more files is listed in both lists."))
  }
} else if (fileOk("reads002.fq") && fileNotOk("reads003.fq")) {
  # Case: no list file, but only two fastq inputs
  in.sorted <- input.names[order(input.names[,2]),]
  reads <- grep("reads",in.sorted[,1],value = TRUE)
  reads1 <- reads[1]
  reads2 <- reads[2]
} else {
  # Case: no list files, more than two fastq inputs
  stop(paste('CHIPSTER-NOTE: ',"List file is missing. You need to provide a list of read files for both directions."))
}

# command
command <- paste(star.binary,"--genomeDir",path.star.index,"--readFilesIn",reads1,reads2,"--outSAMtype BAM SortedByCoordinate","--twopassMode Basic","--runThreadN",chipster.threads.max,"--alignSJoverhangMin 8","--alignSJDBoverhangMin 1","--outSAMstrandField intronMotif","--outFilterType BySJout","--outFilterMultimapNmax",alignments.per.read,"--outFilterMismatchNmax",mismatches.per.pair)
# Use GTF if provided
if (fileOk("annotation.gtf")) {
  command <- paste(command,"--sjdbGTFfile annotation.gtf")
} else {
  gtf.path <- c(file.path(chipster.tools.path,"genomes","gtf"))
  file.list <- list.files(gtf.path,pattern = "\\.gtf$")
  gtf.file <- grep(organism,file.list,value = TRUE)
  gtf.path <- c(file.path(gtf.path,gtf.file[1]))
  command <- paste(command,"--sjdbGTFfile",gtf.path)
}
command <- paste(command,"--alignIntronMin",align.intron.min)
command <- paste(command,"--alignIntronMax",align.intron.max)
command <- paste(command,"--alignMatesGapMax",align.mates.gap.max)
command <- paste(command,"--outFilterMismatchNoverLmax",out.filter.mismatch.nover.lmax)


# Run STAR
documentCommand(command)
runExternal(command)

# rename result files according to the parameter
if (log.files == "final_log") {
  runExternal("mv Log.final.out Log_final.txt")
} else if (log.files == "final_and_progress") {
  runExternal("mv Log.progress.out Log_progress.txt")
  runExternal("mv Log.final.out Log_final.txt")
}
runExternal("mv Aligned.sortedByCoord.out.bam alignment.bam")

# Change file named in BAM header to display names
displayNamesToBAM("alignment.bam", samtools.binary)

# index bam
if (index.file == "index_file") {
  runExternal(paste(samtools.binary,"index alignment.bam"))
}
# Determine base name
inputnames <- read_input_definitions()
basename <- strip_name(inputnames$reads001.fq)

# Make a matrix of output names
outputnames <- matrix(NA,nrow = 2,ncol = 2)
outputnames[1,] <- c("alignment.bam",paste(basename,".bam",sep = ""))
outputnames[2,] <- c("alignment.bam.bai",paste(basename,".bam.bai",sep = ""))

# Write output definitions file
write_output_definitions(outputnames)

# save version information
star.version <- system(paste(star.binary,"--version"),intern = TRUE)
documentVersion("STAR",star.version)

samtools.version <- system(paste(samtools.binary,"--version | grep samtools"),intern = TRUE)
documentVersion("Samtools",samtools.version)
