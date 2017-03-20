# TOOL concatenate.R: "Concatenate files" (Writes the content of selected files into one output file. You can use this tool for example for merging FASTQ files.)
# INPUT input{...}.data: "Input files" TYPE GENERIC
# OUTPUT OPTIONAL concatenated.txt
# OUTPUT OPTIONAL concatenated.html
# OUTPUT OPTIONAL concatenated.tsv
# OUTPUT OPTIONAL concatenated.csv
# OUTPUT OPTIONAL concatenated.fasta
# OUTPUT OPTIONAL concatenated.fastq
# OUTPUT OPTIONAL concatenated.xml
# OUTPUT OPTIONAL concatenated.bed
# OUTPUT OPTIONAL concatenated.gff3
# OUTPUT OPTIONAL concatenated.text
# PARAMETER outformat: "Output name extension" TYPE [ auto: "Automatic", txt: ".txt, General text formatted data", tsv: ".tsv, Tab delimited table", html: ".html, HTML formatted data", fasta: ".fasta, FASTA formatted sequence data",  fastq: ".fastq, Sequence data in FASTQ format " ] DEFAULT auto (Output format type)
# PARAMETER order: "Order by file name" TYPE [yes, no] DEFAULT yes (Concatenate files in alphabetical order by file name.)

#K.M. 18.4. 2016
source(file.path(chipster.common.path, "zip-utils.R"))
source(file.path(chipster.common.path, "tool-utils.R"))

#check if files are zipped
input.names <- read.table("chipster-inputs.tsv", header=F, sep="\t")
for (i in 1:nrow(input.names)) {
	unzipIfGZipFile(input.names[i,1])
}

if (outformat == "auto"){
 library(tools)	
 outformat <- file_ext(input.names[1,2])
 fcheck <- 0
   #check if all files have same name extension
   for (i in 1:nrow(input.names)) {
	 outformat2 <- file_ext(input.names[i,2])
	 if ( outformat2 == "aln" ||outformat == "bai" ||outformat == "bam" ||outformat == "gifasta" ||outformat == "gz" ||outformat == "pdf" ||outformat == "phylip" ||outformat == "png" ||outformat == "qual" || outformat == "svg"){
		 stop("CHIPSTER-NOTE:", outformat2, " Formatted files should not be concatenated with this tool. If you are sure you want to do the concatenation, please select explicit output file name extension.")
	 }
	 
	 if ( outformat != outformat2){
		 fcheck <- 1
	 }		      
   }
 if ( fcheck != 0) {
	 stop("CHIPSTER-NOTE: Input files contain different file name extensions! If you are sure you want to do the concatenation, please select explicit output file name extension.")
 } 
}




if (outformat == "txt" || outformat == "html" || outformat == "tsv" || outformat == "csv" || outformat == "fasta" || outformat == "fastq" || outformat == "xml" || outformat == "bed" || outformat == "gff3" || outformat == "text" ){
 } else { outformat <- "txt" }


# Make list of inputs

input.list <- vector(mode="character", length=0)
 
if (order == "yes"){
    in.sorted <- input.names[order(input.names[,2]),]
    for (i in 1:nrow(input.names)) {
        input.list <- c(input.list, paste(in.sorted[i,1]))
    }
}else{
    for (i in 1:nrow(input.names)) {
	    input.list <- c(input.list, paste(input.names[i,1]))
	}
}
 
files <- paste(input.list, sep="", collapse=" ")

outfile <- paste("concatenated.",outformat , sep="" )
command.full <- paste("cat", files, "> ", outfile)

system(command.full)

