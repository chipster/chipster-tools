# TOOL cufflinks2.R: "Assemble transcripts using Cufflinks" (Given aligned RNA-seq reads as a BAM file, Cufflinks assembles the alignments into a parsimonious set of transcripts. It then estimates the relative abundances of these transcripts based on how many reads support each one. It is recommended to create the input BAM files using the TopHat aligner. You can merge the resulting GTF files of several samples using the Cuffmerge tool, and use the merged GTF file in differential expression analysis using Cuffdiff.)
# INPUT alignment.bam: "BAM file" TYPE BAM
# INPUT OPTIONAL annotation.gtf: "Reference annotation GTF" TYPE GTF
# INPUT OPTIONAL genome.fa: "Genome for bias correction" TYPE GENERIC
# OUTPUT OPTIONAL genes.fpkm_tracking.tsv
# OUTPUT OPTIONAL isoforms.fpkm_tracking.tsv
# OUTPUT OPTIONAL skipped.gtf
# OUTPUT OPTIONAL transcripts.gtf
# PARAMETER chr: "Chromosome names in my BAM file look like" TYPE [chr1, 1] DEFAULT 1 (Chromosome names must match in the BAM file and in the reference annotation. Check your BAM and choose accordingly.)
# PARAMETER OPTIONAL organism: "Reference organism" TYPE [other, "FILES genomes/gtf .gtf"] DEFAULT other (You can use own GTF file or one of those provided on the server.)
# PARAMETER OPTIONAL ug: "Estimate expression of known isoforms, don't assemble novel ones" TYPE [yes: "Only known isoforms", no: "Assemble also novel ones"] DEFAULT no (Reference annotation (a GFF/GTF file\) is used to estimate isoform expression. Program will not assemble novel transcripts, and the it will ignore alignments not structurally compatible with any reference transcript. You can supply your own GTF file or use one of the provided annotations.)
# PARAMETER OPTIONAL lg: "If assembling novel transcripts, should the RABT approach be used?" TYPE [yes, no] DEFAULT no (If you select yes, Cufflinks will do reference annotation based transcript assembly. In RABT a GFF/GTF file is used to guide the assembly, you can supply your own or use one of the provided ones. Reference transcripts will be tiled with faux-reads to provide additional information in assembly. Note that RABT is very slow and it is not recommended for organisms with deep annotations as it can produce false positives. Output will include all reference transcripts as well as any novel genes and isoforms that are assembled.)
# PARAMETER OPTIONAL mmread: "Enable multi-mapped read correction" TYPE [yes, no] DEFAULT no (By default, Cufflinks will uniformly divide each multi-mapped read to all of the positions it maps to. If multi-mapped read correction is enabled, Cufflinks will re-estimate the transcript abundances dividing each multi-mapped read probabilistically based on the initial abundance estimation, the inferred fragment length and fragment bias, if bias correction is enabled.)
# PARAMETER OPTIONAL bias: "Correct for sequence-specific bias" TYPE [yes, no] DEFAULT no (Cufflinks can detect sequence-specific bias and correct for it in abundance estimation. You will need to supply a reference genome as a FASTA file if you are not using one of the provided reference organisms.)


# AMS 21.11.2012
# AMS 2.07.2013 Added chr1/1 option
# EK 3.11.2013 Renamed bias correction parameter, removed quartile normalization parameter
# AMS 11.11.2013 Added thread support
# AMS 24.2.2014 Added G option, changed g option, added support for internal GTFs for both
# AMS 04.07.2014 New genome/gtf/index locations & names

# check out if the file is compressed and if so unzip it
source(file.path(chipster.common.lib.path, "zip-utils.R"))
unzipIfGZipFile("annotation.gtf")
unzipIfGZipFile("genome.fa")

source(file.path(chipster.common.lib.path, "tool-utils.R"))

# binary
cufflinks.binary <- c(file.path(chipster.tools.path, "cufflinks2", "cufflinks"))

# options
cufflinks.options <- ""
cufflinks.options <- paste(cufflinks.options, "-p", chipster.threads.max)

# if (normalize == "yes") {
# 	cufflinks.options <- paste(cufflinks.options, "-N")
# }

# If user has provided a GTF, we use it
if (organism == "other") {
    # If user has provided a GTF, we use it
    if (file.exists("annotation.gtf")) {
        annotation.file <- "annotation.gtf"
    } else {
        stop(paste("CHIPSTER-NOTE: ", "You need provide a GTF file if you are not using one of the provided ones."))
    }
} else {
    # If not, we use the internal one.
    internal.gtf <- file.path(chipster.tools.path, "genomes", "gtf", paste(organism, ".gtf", sep = "", collapse = ""))
    # If chromosome names in BAM have chr, we make a temporary copy of gtf with chr names, otherwise we use it as is.
    if (chr == "chr1") {
        source(file.path(chipster.common.lib.path, "gtf-utils.R"))
        addChrToGtf(internal.gtf, "internal_chr.gtf")
        annotation.file <- paste("internal_chr.gtf")
    } else {
        annotation.file <- paste(internal.gtf)
    }
}

if (lg == "yes") {
    cufflinks.options <- paste(cufflinks.options, "-g", annotation.file)
}
if (ug == "yes") {
    cufflinks.options <- paste(cufflinks.options, "-G", annotation.file)
}

if (mmread == "yes") {
    cufflinks.options <- paste(cufflinks.options, "-u")
}

if (bias == "yes") {
    if (organism == "other") {
        # If user has provided a FASTA, we use it
        if (file.exists("genome.fa")) {
            refseq <- paste("genome.fa")
        } else {
            stop(paste("CHIPSTER-NOTE: ", "If you choose to use bias correction, you need to provide a genome FASTA."))
        }
    } else {
        # If not, we use the internal one. Because FASTA name may lack the version number GTF name has, we
        # first remove the number and then match the base name to folder listing of FASTA files.
        fasta.folder <- file.path(chipster.tools.path, "genomes", "fasta")
        fasta.listing <- list.files(fasta.folder, pattern = "*.fa$")
        organism.base <- remove_extension(organism)
        fasta.name <- grep(organism.base, fasta.listing, value = TRUE)
        internal.fa <- file.path(fasta.folder, fasta.name)
        # internal.fa <- file.path(chipster.tools.path, "genomes", "fasta", paste(organism, ".fa" ,sep="" ,collapse=""))

        # If chromosome names in BAM have chr, we make a temporary copy of fasta with chr names, otherwise we use it as is.
        if (chr == "chr1") {
            source(file.path(chipster.common.lib.path, "seq-utils.R"))
            addChrToFasta(internal.fa, "internal_chr.fa")
            refseq <- paste("internal_chr.fa")
        } else {
            refseq <- paste(internal.fa)
        }
    }
    cufflinks.options <- paste(cufflinks.options, "-b", refseq)
}

# command
command <- paste(cufflinks.binary, cufflinks.options, "-q", "-o tmp", "alignment.bam")

# run
system(command)

# Rename files
source(file.path(chipster.common.lib.path, "gtf-utils.R"))

if (file.exists("tmp/genes.fpkm_tracking") && file.info("tmp/genes.fpkm_tracking")$size > 0) {
    system("mv tmp/genes.fpkm_tracking genes.fpkm_tracking.tsv")
}
if (file.exists("tmp/isoforms.fpkm_tracking") && file.info("tmp/isoforms.fpkm_tracking")$size > 0) {
    system("mv tmp/isoforms.fpkm_tracking isoforms.fpkm_tracking.tsv")
}
# if (file.exists("tmp/skipped.gtf") && file.info("tmp/skipped.gtf")$size > 0) {
# 	sort.gtf("tmp/skipped.gtf", "skipped.gtf")
# }
if (file.exists("tmp/transcripts.gtf") && file.info("tmp/transcripts.gtf")$size > 0) {
    system("mv tmp/transcripts.gtf transcripts.gtf")
}
