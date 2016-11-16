# TOOL htseq-count.R: "Count aligned reads per genes with HTSeq" (Calculates how many reads in a BAM file map to each gene. If you have stranded data, please read the description for strandedness options carefully. You can find more information also in the manual. This tool is based on the HTSeq package. In order to use the output count files for differential expression analysis in edgeR or DESeq, you need to select all the samples and run the tool \"Utilities - Define NGS experiment\".)
# INPUT alignment.bam: "BAM alignment file" TYPE GENERIC
# OUTPUT htseq-counts.tsv
# OUTPUT OPTIONAL htseq-count-info.txt
# PARAMETER organism: "Reference organism" TYPE [Arabidopsis_thaliana.TAIR10.32, Bos_taurus.UMD3.1.85, Canis_familiaris.BROADD2.67, Canis_familiaris.CanFam3.1.85, Drosophila_melanogaster.BDGP5.78, Drosophila_melanogaster.BDGP6.85, Felis_catus.Felis_catus_6.2.85, Gallus_gallus.Galgal4.85, Gasterosteus_aculeatus.BROADS1.85, Halorubrum_lacusprofundi_atcc_49239.ASM2220v1.32, Homo_sapiens.GRCh37.75, Homo_sapiens.GRCh38.85, Homo_sapiens.NCBI36.54, Medicago_truncatula.MedtrA17_4.0.32, Mus_musculus.GRCm38.85, Mus_musculus.NCBIM37.67, Oryza_sativa.IRGSP-1.0.32, Ovis_aries.Oar_v3.1.85, Populus_trichocarpa.JGI2.0.32, Rattus_norvegicus.RGSC3.4.69, Rattus_norvegicus.Rnor_5.0.79, Rattus_norvegicus.Rnor_6.0.85, Schizosaccharomyces_pombe.ASM294v2.32, Solanum_tuberosum.SolTub_3.0.32, Sus_scrofa.Sscrofa10.2.85, Vitis_vinifera.IGGP_12x.32, Yersinia_enterocolitica_subsp_palearctica_y11.ASM25317v1.32, Yersinia_pseudotuberculosis_ip_32953_gca_000834295.ASM83429v1.32] DEFAULT Homo_sapiens.GRCh38.85 (Which organism is your data from.)
# PARAMETER chr: "Chromosome names in the BAM file look like" TYPE [chr1, 1] DEFAULT 1 (Chromosome names must match in the BAM file and in the reference annotation. Check your BAM and choose accordingly.)
# PARAMETER paired: "Does the BAM file contain paired-end data" TYPE [yes, no] DEFAULT no (Does the alignment data contain paired end or single end reads?)
# PARAMETER stranded: "Is the data stranded and how" TYPE [reverse:"\"reverse\" in HTSeq\: the second read of a pair should map to the same strand as the gene", yes:"\"yes\" in HTSeq\: the first read should map to the same strand as the gene", no:"\"no\" in HTSeq\: the data is unstranded"] DEFAULT no (If you select NO, a read will be counted for a gene regardless of which strand it maps to. If you select YES and you have single end data, the read has to map to the same strand as the gene. For paired end data, the first read of a pair has to map to the same strand as the gene, and the second read has to map to the opposite strand. If you select REVERSE and you have paired end data, the second read has to map to the same strand as the gene, and the first read has to map to the opposite strand. You should use REVERSE for paired end data produced for example with the Illumina TruSeq Stranded kit.)
# PARAMETER OPTIONAL mode: "Mode to handle reads overlapping more than one feature" TYPE [union, intersection-strict, intersection-nonempty] DEFAULT union (How to deal with reads that overlap more than one gene or exon?)
# PARAMETER OPTIONAL minaqual: "Minimum alignment quality" TYPE INTEGER FROM 0 TO 100 DEFAULT 10 (Skip all reads with alignment quality lower than the given minimum value.)
# PARAMETER OPTIONAL feature.type: "Feature type to count" TYPE [exon, CDS] DEFAULT exon (Which feature type to use, all features of other type are ignored.)
# PARAMETER OPTIONAL id.attribute: "Feature ID to use" TYPE [gene_id, transcript_id, gene_name, transcript_name, protein_name] DEFAULT gene_id (GFF attribute to be used as feature ID. Several GFF lines with the same feature ID will be considered as parts of the same feature. The feature ID is used to identify the counts in the output table.)
# PARAMETER OPTIONAL print.coord: "Add chromosomal coordinates to the count table" TYPE [yes, no] DEFAULT yes (If you select yes, chromosomal coordinates are added to the output file. Given are the minimum and maximum coordinates of features, e.g. exons, associated with a given identifier)

# 18.1.2012 TH and EK 
# 17.4.2012 EK changed to use Ensembl GTFs 
# 3.2.2013 AMS added chr/nochr option
# 6.5.2013 MK added chr-location information to the output
# 30.5.2013 EK changed the default for "add chromosomal coordinates" to no
# 21.5.2014 EK updated to use HTSeq 0.6.1
# 19.6.2014 AMS changed handling of GTFs
# 4.07.2014 AMS New genome/gtf/index locations & names
# 22.9.2016 EK Clarified strandedness options

# sort bam if the data is paired-end
samtools.binary <- file.path(chipster.tools.path, "samtools", "samtools")
if(paired == "yes"){
	system(paste(samtools.binary, "sort -n alignment.bam name-sorted"))
	bam<-"name-sorted.bam"
} else {
	bam<-"alignment.bam"
}

# htseq-count
if(print.coord == "no") {
	htseq.binary <- file.path(chipster.tools.path, "htseq", "htseq-count")
} else {
	htseq.binary <- file.path(chipster.tools.path, "htseq", "htseq-count_chr")
}


internal.gtf <- file.path(chipster.tools.path, "genomes", "gtf", paste(organism, ".gtf" ,sep="" ,collapse=""))
if(chr == "1"){
	annotation.file <- paste(internal.gtf)
}else{
	source(file.path(chipster.common.path, "gtf-utils.R"))
	addChrToGtf(internal.gtf, "internal_chr.gtf") 
	annotation.file <- paste("internal_chr.gtf")
}


htseq <- paste(htseq.binary, "-f bam -q -m", mode, "-s", stranded, "-a", minaqual, "-t", feature.type, "-i", id.attribute, bam, annotation.file, " > htseq-counts-out.txt")

# run
system(htseq)

# separate result file
system("head -n -5 htseq-counts-out.txt > htseq-counts.tsv")
system("tail -n 5 htseq-counts-out.txt > htseq-count-info.txt")

# bring in file to R environment for formating
file <- c("htseq-counts.tsv")
dat <- read.table(file, header=F, sep="\t")

if(print.coord == "no") {
	names(dat) <- c("id", "count")
} else {
	names(dat) <- c("id", "chr", "start", "end", "len", "strand", "count")
}

# write result table to output
write.table(dat, file="htseq-counts.tsv", col.names=T, quote=F, sep="\t", row.names=F)

# Handle output names
source(file.path(chipster.common.path, "tool-utils.R"))

# read input names
inputnames <- read_input_definitions()

basename <- strip_name(inputnames$alignment.bam)

# Make a matrix of output names
outputnames <- matrix(NA, nrow=1, ncol=2)
outputnames[1,] <- c("htseq-counts.tsv", paste(basename, ".tsv", sep =""))

# Write output definitions file
write_output_definitions(outputnames)

# EOF


