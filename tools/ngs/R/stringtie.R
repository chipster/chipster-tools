# TOOL stringtie.R: "Assemble transcripts using Stringtie" (Stringtie)
# INPUT alignment.bam: "BAM file" TYPE BAM
# OUTPUT OPTIONAL out.gtf

# INPUT OPTIONAL annotation.gtf: "Reference annotation GTF" TYPE GTF
# INPUT OPTIONAL genome.fa: "Genome for bias correction" TYPE GENERIC
# OUTPUT OPTIONAL genes.fpkm_tracking.tsv  
# OUTPUT OPTIONAL isoforms.fpkm_tracking.tsv
# OUTPUT OPTIONAL skipped.gtf
# OUTPUT OPTIONAL transcripts.gtf
# PARAMETER chr: "Chromosome names in my BAM file look like" TYPE [chr1, 1] DEFAULT 1 (Chromosome names must match in the BAM file and in the reference annotation. Check your BAM and choose accordingly.)
# PARAMETER OPTIONAL organism: "Reference organism" TYPE [other, "FILES genomes/gtf .gtf"] DEFAULT other (You can use own GTF file or one of those provided on the server.)
# PARAMETER OPTIONAL ug: "Estimate expression of known isoforms, don't assemble novel ones" TYPE [yes, no] DEFAULT no (Reference annotation (a GFF/GTF file\) is used to estimate isoform expression. Program will not assemble novel transcripts, and the it will ignore alignments not structurally compatible with any reference transcript. You can supply your own GTF file or use one of the provided annotations.)
# PARAMETER OPTIONAL lg: "Do reference annotation based transcript assembly" TYPE [yes, no] DEFAULT no (Cufflinks will use the supplied reference annotation (a GFF/GTF file\) to guide RABT assembly. Reference transcripts will be tiled with faux-reads to provide additional information in assembly. Output will include all reference transcripts as well as any novel genes and isoforms that are assembled. You can supply your own GTF file or use one of the provided annotations.)
# PARAMETER OPTIONAL mmread: "Enable multi-mapped read correction" TYPE [yes, no] DEFAULT no (By default, Cufflinks will uniformly divide each multi-mapped read to all of the positions it maps to. If multi-mapped read correction is enabled, Cufflinks will re-estimate the transcript abundances dividing each multi-mapped read probabilistically based on the initial abundance estimation, the inferred fragment length and fragment bias, if bias correction is enabled.)
# PARAMETER OPTIONAL bias: "Correct for sequence-specific bias" TYPE [yes, no] DEFAULT no (Cufflinks can detect sequence-specific bias and correct for it in abundance estimation. You will need to supply a reference genome as a FASTA file if you are not using one of the provided reference organisms.)

# AO 31.7.2017 Init

# check out if the file is compressed and if so unzip it
source(file.path(chipster.common.path, "zip-utils.R"))
unzipIfGZipFile("annotation.gtf")
unzipIfGZipFile("genome.fa")

source(file.path(chipster.common.path, "tool-utils.R"))

# binary
stringtie.binary <- c(file.path(chipster.tools.path, "stringtie", "stringtie"))

# options
stringtie.options <- ""
stringtie.options <- paste(stringtie.options, "-p", chipster.threads.max)

# command
command <- paste(stringtie.binary, "alignment.bam", stringtie.options, "-o out.gtf")

# run
system(command)


