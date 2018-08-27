# TOOL md5sum.R: "Calculate an md5sums " (Calculates md5sum for a file. Alternativevely if two files are given: file and corresponding md5 sum file, checks that the file matches the MD5 checksum)
# INPUT file: "File" TYPE GENERIC
# INPUT OPTIONAL md5file: "MD5 file" TYPE GENERIC
# OUTPUT OPTIONAL file.md5

source(file.path(chipster.common.path, "tool-utils.R"))

# We need to change the file name before compressing so it is preserved when uncompressing
input.names <- read.table("chipster-inputs.tsv", header=F, sep="\t")


if(file.exists("md5file")){
    #rename input file
	ls <- system ("ls -l >> md5sum.log" )
	md5_a <- tools::md5sum("file")
	md5_as <- toString(md5_a)
	con = file("md5file", "r")
	a <- readLines(con, n = 1)
	fields = strsplit(a, " ")
	md5_b <- fields[[1]][[1]]
	if( md5_as == md5_b  ){
		stop("CHIPSTER-NOTE: OK. The file matches to the given MD5 sum")
	} else{
		in1 <- input.names[1,2]
		in2 <- input.names[2,2]
		stop(paste ("CHIPSTER-NOTE: ERROR: The file", in1, "  does not match the MD5 sum in file",in2 ) )
	}	
}else{
	# Gzip
	system(paste("mv", input.names[1,1], input.names[1,2]))
    system(paste("md5sum", input.names[1,2], "> file.md5"))

    # 
    filename <- paste(input.names[1,2], ".md5", sep = "")
   

    # And finally change Chipster display name to match original file name
    outputnames <- matrix(NA, nrow=1, ncol=2)
    outputnames[1,] <- c("file.md5", filename)

    # Write output definitions file
    write_output_definitions(outputnames)
}


