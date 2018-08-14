# TOOL fixed_text_conversions.R: "Convert contig names to Mothur format" (If you have combined paired MiSeq reads to contigs in another software, the contig names might not work in Mothur. Given a contig fasta file, this tool converts the contig names to Mothur format.)
# INPUT input.dat: input.dat TYPE GENERIC 
# OUTPUT OPTIONAL output.dat 
# OUTPUT OPTIONAL output.txt 
# OUTPUT OPTIONAL output.fasta 
# PARAMETER OPTIONAL type: "Conversion type" TYPE [clc-mothur: "CLC to Mothur"] DEFAULT clc-mothur (File coversion to be executed)

# check out if the file is compressed and if so unzip it


# Handle output names
source(file.path(chipster.common.path, "tool-utils.R"))
# read input names
inputnames <- read_input_definitions()




source(file.path(chipster.common.path, "zip-utils.R"))
input.names <- read.table("chipster-inputs.tsv", header=F, sep="\t")
for (i in 1:nrow(input.names)) {
	unzipIfGZipFile(input.names[i,1])	
}

if ( type == "clc-mothur" ){
	emboss.path <- file.path(chipster.tools.path, "emboss" ,"bin")	
    #check sequece file type
	inputfile.to.check <- ("query.fa")
#	sfcheck.binary <- file.path(chipster.module.path ,"/shell/sfcheck.sh")
#	Add /../misc to path so that sfcheck.sh is found also when this tool runs in NGS module
	sfcheck.binary <- file.path(chipster.module.path ,"/../misc/shell/sfcheck.sh")
	sfcheck.command <- paste(sfcheck.binary, emboss.path, inputfile.to.check )
	str.filetype <- system(sfcheck.command, intern = TRUE )
	
	if ( str.filetype == "Not an EMBOSS compatible sequence file"){
		stop("CHIPSTER-NOTE: Your input file is not a sequence file that is compatible with the tool you try to use")
	}

	system (" awk -F '_1:N:0' '{ print $1 }' input.dat | sed -e s/:/_/g > output.fasta ") 	
}

# Handle output names
#
# read input names
inputnames <- read_input_definitions()

# Determine base name

# Default name if neither special case match. 
# Special cases are for compatibility with older versions of the script
basename <- strip_name(inputnames$input.dat)


# Make a matrix of output names
outputnames <- matrix(NA, nrow=1, ncol=2)
outputnames[1,] <- c("output.fasta", paste(basename, "_converted.fasta", sep =""))

# Write output definitions file
write_output_definitions(outputnames)

