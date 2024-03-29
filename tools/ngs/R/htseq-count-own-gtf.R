# TOOL htseq-count-own-gtf.R: "Count aligned reads per genes with HTSeq using own GTF" (Calculates how many reads in a BAM file map to each gene. You have to provide the gene locations in the GTF format. Please note that the chromosome names have to be same in the GTF and BAM files. If you have stranded data, please read the description for strandedness options carefully. You can find more information also in the manual. This tool is based on the HTSeq package. In order to use the output count files for differential expression analysis in edgeR or DESeq, you need to select all the samples and run the tool \"Utilities - Define NGS experiment\".)
# INPUT alignment.bam: "BAM alignment file" TYPE GENERIC
# INPUT features.gtf: "GTF feature file" TYPE GENERIC
# OUTPUT htseq-counts.tsv
# OUTPUT OPTIONAL htseq-count-info.txt
# PARAMETER paired: "Does the BAM file contain paired-end data" TYPE [yes, no] DEFAULT no (Does the alignment data contain paired end or single end reads?)
# PARAMETER stranded: "Is the data stranded and how" TYPE [reverse:"\"reverse\" in HTSeq\: the second read of a pair should map to the same strand as the gene", yes:"\"yes\" in HTSeq\: the first read should map to the same strand as the gene", no:"\"no\" in HTSeq\: the data is unstranded"] DEFAULT no (If you select NO, a read will be counted for a gene regardless of which strand it maps to. If you select YES and you have single end data, the read has to map to the same strand as the gene. For paired end data, the first read of a pair has to map to the same strand as the gene, and the second read has to map to the opposite strand. If you select REVERSE and you have paired end data, the second read has to map to the same strand as the gene, and the first read has to map to the opposite strand. You should use REVERSE for paired end data produced for example with the Illumina TruSeq Stranded kit.)
# PARAMETER OPTIONAL mode: "Mode to handle reads overlapping more than one gene" TYPE [union, intersection-strict, intersection-nonempty] DEFAULT union (How to deal with reads that overlap more than one gene or exon?)
# PARAMETER OPTIONAL minaqual: "Minimum alignment quality" TYPE INTEGER FROM 0 TO 100 DEFAULT 10 (Skip all reads with alignment quality lower than the given minimum value.)
# PARAMETER OPTIONAL feature.type: "Feature type to count" TYPE [exon, CDS, gene] DEFAULT exon (Which feature type to use, all features of other type are ignored.)
# PARAMETER OPTIONAL id.attribute: "Feature ID to use" TYPE [gene_id, ID, GeneID, transcript_id, gene_name, transcript_name, protein_name, proteinId] DEFAULT gene_id (GFF attribute to be used as feature ID. Several GFF lines with the same feature ID will be considered as parts of the same feature. The feature ID is used to identity the counts in the output table.)
# PARAMETER OPTIONAL print.coord: "Add chromosomal coordinates to the count table" TYPE [yes, no] DEFAULT yes (If you select yes, chromosomal coordinates are added to the output file. Given are the minimum and maximum coordinates of features, e.g. exons, associated with a given identifier)


# 22.8.2011	TH and EK
# 6.5.2013	MK added chr-location information to the output
# 21.5.2014	EK updated to use HTSeq 0.6.1
# 9.4.2015	ML added the geneID options
# 22.9.2016	EK clarified strandedness options
# 24.3.2021	EK added feature type gene
# 24.3.2021	EK added id.attribute proteinId

# check out if the file is compressed and if so unzip it
source(file.path(chipster.common.lib.path, "zip-utils.R"))
unzipIfGZipFile("features.gtf")

# sort bam if the data is paired-end
samtools.binary <- file.path(chipster.tools.path, "samtools-0.1.19", "samtools")
if (paired == "yes") {
    system(paste(samtools.binary, "sort -n alignment.bam name-sorted"))
    bam <- "name-sorted.bam"
} else {
    bam <- "alignment.bam"
}

# htseq-count
if (print.coord == "no") {
    htseq.binary <- file.path(chipster.tools.path, "htseq", "htseq-count")
} else {
    htseq.binary <- file.path(chipster.tools.path, "htseq", "htseq-count_chr")
}

htseq <- paste(htseq.binary, "-f bam -q -m", mode, "-s", stranded, "-a", minaqual, "-t", feature.type, "-i", id.attribute, bam, "features.gtf > htseq-counts-out.txt")

# run
system(htseq)

# separate result file
system("head -n -5 htseq-counts-out.txt > htseq-counts.tsv")
system("tail -n 5 htseq-counts-out.txt > htseq-count-info.txt")

# bring in file to R environment for formating
file <- c("htseq-counts.tsv")
dat <- read.table(file, header = F, sep = "\t")
if (print.coord == "no") {
    names(dat) <- c("id", "count")
} else {
    names(dat) <- c("id", "chr", "start", "end", "len", "strand", "count")
}

# write result table to output
write.table(dat, file = "htseq-counts.tsv", col.names = T, quote = F, sep = "\t", row.names = F)

# Add additional info lines about read totals to output
file2 <- c("htseq-count-info.txt")
dat2 <- read.table(file2, header = F, sep = "\t")

assigned <- sum(dat$count)
notassigned <- sum(dat2[2])
total <- assigned + notassigned

line <- paste("\n")
line <- paste(line, "not_counted\t", notassigned, "\n", sep = "")
line <- paste(line, "counted\t", assigned, "\n", sep = "")
line <- paste(line, "total\t", total, "\n", sep = "")

write(line, "htseq-count-info.txt", append = TRUE)

# Handle output names
source(file.path(chipster.common.lib.path, "tool-utils.R"))

# read input names
inputnames <- read_input_definitions()

basename <- strip_name(inputnames$alignment.bam)

# Make a matrix of output names
outputnames <- matrix(NA, nrow = 1, ncol = 2)
outputnames[1, ] <- c("htseq-counts.tsv", paste(basename, ".tsv", sep = ""))

# Write output definitions file
write_output_definitions(outputnames)

# EOF
