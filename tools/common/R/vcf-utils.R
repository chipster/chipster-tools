# Check VCF chromosome naming scheme
getVCFNames <- function(input){
    # Uncompress
	unzipIfGZipFile(input)

    # Subset file to make searching faster
    system(paste("grep -v \"#\"", input, "|head -50 > name.tmp"))
    
    # Check if VCF has chr chromosome names
    vcf.chr <- grep("^chr", readLines("name.tmp"), ignore.case = TRUE)
    if (length(vcf.chr) > 0){
        vcf.names <- paste("chr1")
    }else{
        vcf.names <- paste("1")
    }
    return(vcf.names)
}

# Adds "chr" to the beginning of each line of a VCF file that starts with a number or with X, Y, Z, W or M
#
addChrToVCF <- function(input, output){
    system(paste("sed /^.[0-9,X,Y,Z,W,M]/s/\">\"/\">chr\"/  <", input, ">", output))
}

# Removes "chr"  from each of a VCF file line that starts wirh "chr" (case insensitve).
#
removeChrFromVCF <- function(input, output){
    #system(paste("awk '{gsub(/^chr/,""); print}'", input, ">", output
    system(paste("cat", input, "|sed s/ID=chr/ID=/I | sed s/^chr//I >", output))
}
 