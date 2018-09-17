# TOOL muscle_profile.R: "Add a sequences to an existing alignmnet with MUSCLE" (Add a sequence or multiple sequence alignment to an existing alignment with MUSCLE software.)
# INPUT sequence1: "Alignmnet" TYPE GENERIC
# INPUT sequence2: "Sequence of alignmnet to add" TYPE GENERIC
# OUTPUT OPTIONAL alignment.aln.txt
# OUTPUT OPTIONAL alignment.html 
# OUTPUT OPTIONAL alignment.fasta 
# OUTPUT OPTIONAL alignment.log
# PARAMETER OPTIONAL maxiters: "Maximum iterations" TYPE INTEGER DEFAULT 16 (Maximum number of iterations )
# PARAMETER OPTIONAL outformat: "Output format" TYPE [ aln: "Clustal", fasta: "Fasta", html: "HTML" ] DEFAULT fasta (Multiple sequence alingment file format)
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


## These are for permanent conda-muscle installation
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
seqcount.exe <- file.path(emboss.path, "seqcount -filter sequence2")
str.queryseq <- system(seqcount.exe, intern = TRUE )
num.queryseq <- as.integer(str.queryseq)
#round(num.queryseq)

if (num.queryseq > 100000){
	stop(paste('CHIPSTER-NOTE: Too many query sequences. Maximun is 100000 but your file contains ', num.queryseq ))
}


if ( str.filetype != "fasta"){
	seqret.binary <- file.path(emboss.path,"seqret")
	seqret.command <- paste(seqret.binary, "sequence2 sequence2.fasta -auto")
	system(seqret.command)
	system("rm -f sequence2")
	system("mv sequence2.fasta sequence2")
}

#count the query sequeces
seqcount.exe <- file.path(emboss.path, "seqcount -filter sequence1")
str.queryseq <- system(seqcount.exe, intern = TRUE )
num.queryseq <- as.integer(str.queryseq)
#round(num.queryseq)

if (num.queryseq > 100000){
	stop(paste('CHIPSTER-NOTE: Too many query sequences. Maximun is 100000 but your file contains ', num.queryseq ))
}


if ( str.filetype != "fasta"){
	seqret.binary <- file.path(emboss.path,"seqret")
	seqret.command <- paste(seqret.binary, "sequence1 sequence1.fasta -auto")
	system(seqret.command)
	system("rm -f sequence1")
	system("mv sequence1.fasta sequence1")
}

#muscle.binary <- paste(conda.path, "conda_muscle")

muscle.options <- paste("-maxiters", maxiters)


outfile <- ("alignment.fasta")

if (outformat == "aln" ){
	muscle.options <- paste(muscle.options, "-clw")
	outfile <- ("alignment.aln.txt")
}

if (outformat == "html" ){
	muscle.options <- paste(muscle.options, "-html")
	outfile <- ("alignment.html")
}



command.full <- paste(conda.path, conda.def, " -profile -in1 sequence1 -in2 sequence2 ", muscle.options, " -out", outfile, "2>alignment.log" )
system(command.full)

if ( save_log == "no") {
	system ("rm -f alignment.log")
}
