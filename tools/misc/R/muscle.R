# TOOL muscle.R: "Multiple sequence alignment with MUSCLE" (Multiple sequence alignmnet with MUSCLE software.)
# INPUT sequence: "Query sequences" TYPE GENERIC 
# OUTPUT OPTIONAL alignment.aln.txt
# OUTPUT OPTIONAL alignment.html 
# OUTPUT OPTIONAL alignment.fasta 
# OUTPUT OPTIONAL newick_tree.txt
# OUTPUT OPTIONAL alignment.log
# PARAMETER OPTIONAL maxiters: "Maximum iterations" TYPE INTEGER DEFAULT 16 (Maximum number of iterations )
# PARAMETER OPTIONAL outformat: "Output format" TYPE [ aln: "Clustal", fasta: "Fasta"] DEFAULT fasta (Multiple sequence alingment file format)
# PARAMETER OPTIONAL html: "HTML formatted alignmnet" TYPE [ yes: "yes", no: "no" ] DEFAULT no (Print HTML formatted alignment file)
# PARAMETER OPTIONAL njtree: "Print tree" TYPE [ yes: "yes", no: "no" ] DEFAULT no (Create Neighbor-joining thee based on the alignmnet. Treefile is in Newick format)
# PARAMETER OPTIONAL save_log: "Collect a log file" TYPE [yes: yes, no: no] DEFAULT no (Collect a log file.)


source(file.path(chipster.common.path, "zip-utils.R"))
unzipIfGZipFile("sequence")
emboss.path <- file.path(chipster.tools.path, "emboss" ,"bin")

conda_otf.binary <- file.path(chipster.module.path, "../admin/shell/conda_otf.sh ")
#what conda package to install
conda.package <- ("muscle")
#what command to launch
conda.tool <- ("muscle")
conda.def <- paste(conda.package, "/", conda.tool, sep="")

#replacement for a single binary
muscle.binary <- paste(conda_otf.binary, conda.def)

#conda.path <- file.path("/opt/chipster/tools_local" ,"miniconda2","conda_execute")
#conda.env <- ("chipster-phy")
#conda.tool <- ("muscle")
#conda.def <- paste(conda.env, "/", conda.tool, sep="")
#check sequece file type16
sfcheck.binary <- file.path(chipster.module.path ,"/shell/sfcheck.sh")
sfcheck.command <- paste(sfcheck.binary, emboss.path, "sequence" )
str.filetype <- system(sfcheck.command, intern = TRUE )

if ( str.filetype == "Not an EMBOSS compatible sequence file"){
	stop("CHIPSTER-NOTE: Your input file is not a sequence file that is compatible with the tool you try to use")
}

#count the query sequeces
seqcount.exe <- file.path(emboss.path, "seqcount -filter sequence")
str.queryseq <- system(seqcount.exe, intern = TRUE )
num.queryseq <- as.integer(str.queryseq)
#round(num.queryseq)

if (num.queryseq > 100000){
	stop(paste('CHIPSTER-NOTE: Too many query sequences. Maximun is 100000 but your file contains ', num.queryseq ))
}


if ( str.filetype != "fasta"){
	seqret.binary <- file.path(emboss.path,"seqret")
	seqret.command <- paste(seqret.binary, "sequence sequence.fasta -auto")
	system(seqret.command)
	system("rm -f sequence")
	system("mv sequence.fasta sequence")
}


#muscle.binary <- paste(conda.path, "conda_muscle")

muscle.options <- paste("-maxiters", maxiters)


if (outformat == "fasta" ){
	muscle.options <- paste(muscle.options, "-fastaout alignment.fasta ")
	outfile <- ("alignment.fasta")
}

if (outformat == "aln" ){
	muscle.options <- paste(muscle.options, "-clw alignment.aln.txt")
    outfile <- ("alignment.aln.txt")
}

if (html== "yes" ){
	muscle.options <- paste(muscle.options, "-htmlout alignment.html")
}

command.full <- paste(muscle.binary, " -in sequence ", muscle.options, " 2>>alignment.log" )
#command.full <- paste(conda.path, conda.def, " -in sequence ", muscle.options, " 2>>alignment.log" )
system(command.full)

if (njtree == "yes"){
	command.full <- paste(muscle.binary, "-maketree -in", outfile, " -out newick_tree.txt -cluster neighborjoining 2>>alignment.log" )
	system(command.full)
	
}


if ( save_log == "no") {
	system ("rm -f alignment.log")
}
