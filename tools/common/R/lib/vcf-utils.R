# Check VCF chromosome naming scheme
getVCFNames <- function(input) {
    # Uncompress
    unzipIfGZipFile(input)

    # Subset file to make searching faster
    system(paste("grep -v \"#\"", input, "|head -50 > name.tmp"))

    # Check if VCF has chr chromosome names
    vcf.chr <- grep("^chr", readLines("name.tmp"), ignore.case = TRUE)
    if (length(vcf.chr) > 0) {
        vcf.names <- paste("chr1")
    } else {
        vcf.names <- paste("1")
    }
    return(vcf.names)
}

# Adds "chr" to the chromosome names of a VCF file that start with a number or with X, Y, Z, W or M:
# the CHROM column of the records and the IDs of the ##contig header lines, which GATK checks too.
#
addChrToVCF <- function(input, output) {
    system(paste(
        "sed -e 's/^\\([0-9XYZWM]\\)/chr\\1/' -e 's/^\\(##contig=<ID=\\)\\([0-9XYZWM]\\)/\\1chr\\2/'",
        input,
        ">",
        output
    ))
}

# Removes "chr"  from each of a VCF file line that starts wirh "chr" (case insensitve).
#
removeChrFromVCF <- function(input, output) {
    # system(paste("awk '{gsub(/^chr/,""); print}'", input, ">", output
    system(paste("cat", input, "|sed s/ID=chr/ID=/I | sed s/^chr//I >", output))
}
