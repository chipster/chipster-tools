# TOOL samtools-snp-indel-single.R: "Call SNPs and INDELs with SAMtools" (Calls SNPs and short INDELs for one or more diploid individuals. Alignment records are grouped by sample identifiers in @RG header lines. If sample identifiers are absent, each input file is regarded as one sample. You can provide your own reference sequence in FASTA format or choose one of the provided reference genomes. This tool is based on the SAMtools package.)
# INPUT alignment{...}.bam: "Sorted BAM files" TYPE BAM
# INPUT OPTIONAL ownref.fa: "Reference sequence FASTA" TYPE GENERIC
# OUTPUT variants.vcf
# PARAMETER organism: "Reference sequence" TYPE [other, "FILES genomes/fasta .fa"] DEFAULT other (Reference sequence.)
# PARAMETER chr: "Chromosome names in my BAM file look like" TYPE [chr1, 1] DEFAULT 1 (Chromosome names must match in the BAM file and in the reference sequence. Check your BAM and choose accordingly. This only applies to provided reference genomes.)
# PARAMETER ploidy: "Ploidy" TYPE [haploid, diploid] DEFAULT diploid (Ploidy.)
# PARAMETER OPTIONAL mpileup.ui: "Call only SNPs, INDELs not considered" TYPE [yes, no] DEFAULT no (Do not perform INDEL calling.)
# PARAMETER OPTIONAL vcfutils.d: "Minimum read depth" TYPE INTEGER DEFAULT 2 (Minimum read depth.)
# PARAMETER OPTIONAL vcfutils.ud: "Maximum read depth" TYPE INTEGER DEFAULT 1000 (Maximum read depth. Should be adjusted to about twice the average read depth.)
# PARAMETER OPTIONAL mpileup.q: "Minimum mapping quality for an alignment to be used" TYPE INTEGER DEFAULT 0 (Minimum mapping quality for an alignment to be used. Default is 0.)
# PARAMETER OPTIONAL mpileup.uq: "Minimum base quality for a base to be considered" TYPE INTEGER DEFAULT 20 (Minimum base quality for a base to be considered. Default is 20.)
# PARAMETER OPTIONAL mpileup.ub: "Disable probabilistic realignment for computation of BAC" TYPE [yes, no] DEFAULT no (Disable probabilistic realignment for the computation of base alignment quality (BAQ\). BAQ is the Phred-scaled probability of a read base being misaligned. Applying this option greatly helps to reduce false SNPs caused by misalignments.)
# PARAMETER OPTIONAL mpileup.uc: "Downgrading coefficient" TYPE INTEGER DEFAULT 0 (Coefficient for downgrading mapping quality for reads containing excessive mismatches. The recommended value for the older BWA-SW is 50. A zero value disables this functionality.)
# PARAMETER OPTIONAL mpileup.r: "Call variants only for a certain region" TYPE STRING DEFAULT all (Only generate pileup in defined region. Region given as chromosome:start-end, e.g. 20:131505-131550.)
# PARAMETER OPTIONAL mpileup.ud: "Output per sample read depth" TYPE [yes, no] DEFAULT no (Output per-sample read depth, DP.)
# PARAMETER OPTIONAL mpileup.dv: "Output per sample number of high-quality non-reference bases" TYPE [yes, no] DEFAULT no (Output per sample number of high-quality non-reference bases, DV.)
# PARAMETER OPTIONAL vcftools.info.all: "Output all INFO fields" TYPE [yes, no] DEFAULT no (By default only AC, DP, DP4, INDEL and IDV fields are included in the output. If you select yes, all the INFO fields that are produced by default by Samtools are included. These are AC, AN, BQB, DP, DP4, HDB, ICB, IDV, IMF, INDEL, MQ, MQB, MQOF, MQSB, RPB, SGB and VDB.)

# AMS 30.5.2012
# AMS 17.08.2012 Additional parameters as suggested by Jarno
# AMS 24.6.2014 Changed handling of fasta file
# AMS 04.07.2014 New genome/gtf/index locations & names
# EK  02.10.2015 Modified output options.

# To be added later PARAMETER OPTIONAL mpileup.us: "Output per sample strand bias P-value" TYPE [yes, no] DEFAULT no (Output per-sample Phred-scaled strand bias P-value.)

# check out if the file is compressed and if so unzip it
source(file.path(chipster.common.lib.path, "zip-utils.R"))
unzipIfGZipFile("ownref.fa")

# binaries
samtools.binary <- c(file.path(chipster.tools.path, "samtools-1.2", "samtools"))
bcftools.binary <- c(file.path(chipster.tools.path, "bcftools", "bcftools"))
vcfutils.binary <- c(file.path(chipster.tools.path, "bcftools", "vcfutils.pl"))
vcftools.binary <- c(file.path(chipster.tools.path, "vcftools", "bin", "vcftools"))

# If user provided fasta we use it, else use internal fasta
if (organism == "other") {
    # If user has provided a FASTA, we use it
    if (file.exists("ownref.fa")) {
        refseq <- paste("ownref.fa")
    } else {
        stop(paste("CHIPSTER-NOTE: ", "You need to provide a FASTA file or choose one of the provided reference genomes."))
    }
} else {
    # If not, we use the internal one.
    internal.fa <- file.path(chipster.tools.path, "genomes", "fasta", paste(organism, ".fa", sep = "", collapse = ""))
    # If chromosome names in BAM have chr, we make a temporary copy of fasta with chr names, otherwise we use it as is.
    if (chr == "chr1") {
        source(file.path(chipster.common.lib.path, "seq-utils.R"))
        addChrToFasta(internal.fa, "internal_chr.fa")
        refseq <- paste("internal_chr.fa")
    } else {
        refseq <- paste(internal.fa)
    }
}

# mpileup otions
mpileup.options <- paste("mpileup -u")
if (mpileup.r != "all") {
    mpileup.options <- paste(mpileup.options, "-r", mpileup.r)
}
if (mpileup.ub == "yes") {
    mpileup.options <- paste(mpileup.options, "-B")
}
mpileup.options <- paste(mpileup.options, "-C", mpileup.uc)
mpileup.options <- paste(mpileup.options, "-q", mpileup.q)
mpileup.options <- paste(mpileup.options, "-Q", mpileup.uq)
if (mpileup.ud == "yes") {
    mpileup.options <- paste(mpileup.options, "-t DP")
}
# if (mpileup.us == "yes") {
# 	mpileup.options <- paste(mpileup.options, "-t SP")
# }
if (mpileup.dv == "yes") {
    mpileup.options <- paste(mpileup.options, "-t DV")
}
if (mpileup.ui == "yes") {
    mpileup.options <- paste(mpileup.options, "-I")
}


# bcftools options
bcftools.options <- paste("call --variants-only --multiallelic-caller --output-type z")

if (ploidy == "haploid") {
    system("for n in $( ls *.bam ); do printf \"${n}\t1\n\" >> samplefile.txt; done")
    bcftools.options <- paste(bcftools.options, "-S samplefile.txt")
}

# vcfutils options
vcfutils.options <- paste("varFilter", "-d", vcfutils.d, "-D", vcfutils.ud)


# commands
command1 <- paste(samtools.binary, mpileup.options, "-f", refseq, "*.bam", "|", bcftools.binary, bcftools.options, "- > var.raw.bcf")
command2 <- paste(bcftools.binary, "view var.raw.bcf |", vcfutils.binary, vcfutils.options, "> variants.raw.vcf")
if (vcftools.info.all == "yes") {
    command3 <- paste(vcftools.binary, "--vcf variants.raw.vcf --out vcftools --recode --recode-INFO-all")
} else {
    command3 <- paste(vcftools.binary, "--vcf variants.raw.vcf --out vcftools --recode --recode-INFO DP --recode-INFO DP4 --recode-INFO IDV --recode-INFO INDEL --recode-INFO AC")
}

# run
# stop(paste('CHIPSTER-NOTE: ', command1))
system(command1)
system(command2)
system(command3)
system("mv vcftools.recode.vcf variants.vcf")

# Change bam names in VCF to original names
input.names <- read.table("chipster-inputs.tsv", header = F, sep = "\t")
for (i in 1:nrow(input.names)) {
    sed.command <- paste("s/", input.names[i, 1], "/", input.names[i, 2], "/", sep = "")
    system(paste("sed -i", sed.command, "variants.vcf"))
}
