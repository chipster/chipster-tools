# TOOL concatenate.R: "Concatenate files" (Writes the content of seletced files into one output file)
# INPUT input{...}.data: "Input files" TYPE GENERIC
# OUTPUT OPTIONAL result.txt
# OUTPUT OPTIONAL result.html
# OUTPUT OPTIONAL result.tsv
# OUTPUT OPTIONAL result.csv
# OUTPUT OPTIONAL result.fasta
# OUTPUT OPTIONAL result.fastq
# OUTPUT OPTIONAL result.xml
# OUTPUT OPTIONAL result.bed
# OUTPUT OPTIONAL result.gff3
# OUTPUT OPTIONAL result.text
# PARAMETER outformat: "Output name extension" TYPE [ auto: "Automatic", txt: ".txt, General text formatted data", tsv: ".tsv, Tab delimited table", html: ".html, HTML formatted data", fasta: ".fasta, FASTA formatted sequence data",  fastq: ".fastq, Sequence data in FASTQ format " ] DEFAULT auto (Output format type)

source(file.path(chipster.common.path, "zip-utils.R"))
#unzipIfGZipFile("input1")
#unzipIfGZipFile("input2")
# Change bam names in VCF to original names
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
		 stop("CHIPSTER-NOTE:", outformat2, " formatted files should not be concatenated with this tool. If you are sure you want to do the concatenation, please select explicit output file name extension.")
	 }
	 
	 if ( outformat != outformat2){
		 fcheck <- 1
	 }		      #stop("CHIPSTER-NOTE:",input.names[1,2], outformat, input.names[2,2])
   }
 if ( fcheck != 0) {
	 stop("CHIPSTER-NOTE:Input files contain different file name extensions! If you are sure you want to do the concatenation, please select explicit output file name extension.")
 } 
}




if (outformat == "txt" || outformat == "html" || outformat == "tsv" || outformat == "csv" || outformat == "fasta" || outformat == "fastq" || outformat == "xml" || outformat == "bed" || outformat == "gff3" || outformat == "text" ){
 } else { outformat <- "txt" }


outfile <- paste("result.",outformat , sep="" )
command.full <- paste("cat *.data > ", outfile)

system(command.full)




