# Formats a VCF file to be comatible with GATK.
# Most GATK tools require VCF files to be bgzip compressed and tabix indexed.
#
# Note that e.g. for input file "example.vcf" you end up with files "example.vcf.gz" and "example.vcf.gz.tbi".
# Remember to use correct filenames in GATK commands.
#
formatGatkVcf <- function(input.vcf, bam.names) {
    # source(file.path(chipster.common.lib.path, "zip-utils.R"))
    tabix.binary <- c(file.path(chipster.tools.path, "tabix", "tabix"))
    bgzip.binary <- c(file.path(chipster.tools.path, "tabix", "bgzip"))

    # Uncompress
    unzipIfGZipFile(input.vcf)

    if (!(missing(bam.names))) {
        # Check if VCF chromosome names match BAM names
        vcf.names <- getVCFNames(input.vcf)
        if ((bam.names == "1") && (vcf.names == "chr1")) {
            # BAM names have "1"- Remove "chr" from VCF names.
            removeChrFromVCF(input.vcf, "output.vcf")
            system(paste("mv output.vcf", input.vcf))
        }
        if ((bam.names == "chr1") && (vcf.names == "1")) {
            # BAM names have "chr". Add "chr" to VCF names.
            addChrToVCF(input.vcf, "output.vcf")
            system(paste("mv output.vcf", input.vcf))
        }
    }

    # GATK requires files to compressed with bgzip and indexed
    # Bgzip
    system(paste(bgzip.binary, input.vcf))
    # Index
    gz.name <- paste(input.vcf, ".gz", sep = "")
    system(paste(tabix.binary, "-p vcf", gz.name))
}

# Formast a FASTA format reference sequence to be compatible with GATK
# Most GATK tools require the reference sequence to be indexed and have a dictionary file
#
# Note that e.g for input file "example.fa" you end up with files "example.fa.fai" and "example.fa.dict".
#
formatGatkFasta <- function(fasta) {
    source(file.path(chipster.common.lib.path, "zip-utils.R"))
    picard.binary <- c(file.path(chipster.tools.path, "picard-tools", "picard.jar"))
    samtools.binary <- c(file.path(chipster.tools.path, "samtools-0.1.19", "samtools"))
    # Uncompress
    unzipIfGZipFile(fasta)
    # Index
    system(paste(samtools.binary, "faidx", fasta))
    # Create dictionary file
    system(paste("java -jar", picard.binary, "CreateSequenceDictionary", paste("R=", fasta, sep = "", collapse = ""), paste("O=", fasta, ".dict", sep = "", collapse = "")))
}
