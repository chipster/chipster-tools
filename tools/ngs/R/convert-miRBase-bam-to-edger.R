# TOOL convert-miRBase-bam-to-edger.R: "Convert miRBased BAM file to count table" (This tool takes a BAM file as an input, calculates the number of times each miRNA is identified, and removes the ones for which the count is under the user defined threshold.)
# INPUT bam_file.bam: "Alignment against miRBase in BAM format" TYPE GENERIC
# OUTPUT miRNA-counts.tsv: "A count file suitable for differential expression analysis"
# PARAMETER count_limit: "Count limit" TYPE INTEGER FROM 0 TO 1000 DEFAULT 10 (Keep miRNAs which have more reads than this.)

# EK 07.07.2011
# MK 13.05.2013, fix bug in header formats
# MK 12.05.2014, added check for file size

source(file.path(chipster.common.path, "tool-utils.R"))

# Convert to SAM, grep miRNA name (removes unaligned reads), count, filter on tag number
samtools.binary <- c(file.path(chipster.tools.path,"samtools","samtools"))
counts.command <- paste(samtools.binary,"view -q 10 bam_file.bam | awk '{print $3}' | grep -v '^*' | sort -k1 | uniq -c | awk '{if($1>",count_limit,")print $2\"\t\"$1}'> counts.tsv")
system(counts.command)

input.file <- "counts.tsv"
if (file.info(input.file)$size > 0) {
  dat <- read.table(input.file,header = FALSE,sep = "\t",row.names = NULL)
  colnames(dat) <- c("id","count")
  write.table(data.frame(dat),file = "miRNA-counts.tsv",col.names = TRUE,quote = FALSE,sep = "\t",row.names = FALSE)
} else {
  stop("CHIPSTER-NOTE: No miRNA in your BAM file had the minimum number of reads required.")
}

# Handle output names
# Read input names and strip file extension
inputnames <- read_input_definitions()
input1name <- inputnames$bam_file.bam # name from the unmapped.bam input
input1namestripped <-strip_name(input1name)
outputnames <- matrix(NA, nrow=1, ncol=2)
outputnames[1,] <- c("miRNA-counts.tsv", paste(input1namestripped, "_counts.tsv", sep = ""))
write_output_definitions(outputnames)
# EOF

