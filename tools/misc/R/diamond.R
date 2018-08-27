# TOOL diamond.R: "DIMAOND protein sequence similarity search" ( Sequence similarity search with DIAMONDsoftware.)
# INPUT sequence: "Query sequences" TYPE GENERIC 
# INPUT OPTIONAL reference: "Reference proteins" TYPE GENERIC 
# OUTPUT OPTIONAL diamond.txt
# OUTPUT OPTIONAL diamond.sam 
# OUTPUT OPTIONAL diamond.csv
# OUTPUT OPTIONAL diamond.log
# OUTPUT OPTIONAL unaligned.fasta
# OUTPUT OPTIONAL reference.dmnd
# PARAMETER OPTIONAL db: "Reference database" TYPE [ swiss: "SwissProt", trembl: "TrEMBL", nr: "NR"] DEFAULT swiss (Reference database. This parameter is ignored if referense sequece set is give as file)
# PARAMETER OPTIONAL maxtargetseqs: "Maximum number of hits" TYPE INTEGER DEFAULT 500 (Maximum number hits to report for one query sequence)
# PARAMETER OPTIONAL evalue: "e-value" TYPE DECIMAL DEFAULT 0.001 (E-value threshold )
# PARAMETER OPTIONAL outformat: "Output format" TYPE [ 0: "BLAST pairwise", 5: "BLAST XML", 6: "BLAST tabular", 100: "DIAMOND alignment archive", 101: "SAM"] DEFAULT 0 (Multiple sequence alingment file format)
# PARAMETER OPTIONAL unal: "Report unaligned sequences" TYPE [ yes: "yes", no: "no" ] DEFAULT no (Collect to a sequences file those query sequences that did not had any matches)
# PARAMETER OPTIONAL smode: "Search mode" TYPE [ fast: "Fast", sensitive: "Sensitive", moresensitive: "More sensitive" ] DEFAULT fast (Search mode)
# PARAMETER OPTIONAL matrix: "Matrix" TYPE [BLOSUM45: "BLOSUM45", BLOSUM50: "BLOSUM50", BLOSUM62: "BLOSUM62", BLOSUM80: "BLOSUM80", BLOSUM90: "BLOSUM90"] DEFAULT BLOSUM62 (Weight matrix assigns a score for aligning pairs of residues, and determines overall alignment score. Experimentation has shown that the BLOSUM62 matrix is among the best for detecting most weak protein similarities. For particularly long and weak alignments, the BLOSUM45 matrix may prove superior. For proteins, shorter than 85 residues, the BLOSUM80 matrix may provide better hits"  )
# PARAMETER OPTIONAL save_log: "Collect a log file" TYPE [yes: yes, no: no] DEFAULT no (Collect a log file.)
# PARAMETER OPTIONAL keep_index: "Save the DIAMOND index" TYPE [yes: yes, no: no] DEFAULT no (Write the DIAMOND indexes of the refrence paritein set to a file for re-use.)

use_remote_index <- ("yes")
keep_index <- ("no")
index_url <- ("https://diamond.object.pouta.csc.fi")
source(file.path(chipster.common.path, "zip-utils.R"))
unzipIfGZipFile("sequence")
emboss.path <- file.path(chipster.tools.path, "emboss" ,"bin")

conda_otf.binary <- file.path(chipster.module.path, "../admin/shell/conda_otf.sh ")
#what conda package to install
conda.package <- ("diamond")
#what command to launch
conda.tool <- ("diamond")
conda.def <- paste(conda.package, "/", conda.tool, sep="")

#replacement for a single binary
diamond.binary <- paste(conda_otf.binary, conda.def)

#conda.path <- file.path("/opt/chipster/tools_local" ,"miniconda2","conda_execute")
#conda.env <- ("chipster-phy")
#conda.tool <- ("muscle")
#conda.def <- paste(conda.env, "/", conda.tool, sep="")
#check sequece file type16
sfcheck.binary <- file.path(chipster.module.path ,"/shell/sfcheck.sh")
sfcheck.command <- paste(sfcheck.binary, emboss.path, "sequence" )
str.filetype <- system(sfcheck.command, intern = TRUE )

if ( str.filetype == "Not an EMBOSS compatible sequence file"){
	stop("CHIPSTER-NOTE: Your query sequence file is not a sequence file that is compatible with the tool you try to use")
}


if(file.exists("reference")){
    system("echo Using user provided reference file >> diamond.log")
	use_remote_index <- ("no")
	db <- ("reference")
	system ("mv reference reference.dmnd")
	dbinfo.command <- paste( diamond.binary, "dbinfo --db reference.dmnd | grep -c version")
	check.number <- system(dbinfo.command, intern = TRUE )
	if ( check.number == "1" ){
		system("echo Referense sequence set is already indexed >> diamond.log")
	}    
	if ( check.number != "1" ){
	   system("echo Calculating Diamond indekses >> diamond.log")
	   sfcheck.command <- paste(sfcheck.binary, emboss.path, "reference.dmnd" )
	   str.filetype <- system(sfcheck.command, intern = TRUE )
	
       if ( str.filetype == "Not an EMBOSS compatible sequence file"){
	    	stop("CHIPSTER-NOTE: Your sequence file is not a sequence file that is compatible with the tool you try to use")
	   }
	   system("mv reference.dmnd reference.fasta ")
	   #Check if the query is protein or nucleotode
	   infoseq.exe <- file.path(emboss.path, "infoseq -nowarning -noerror -nodie -only -type -nohead -filter | uniq | awk '{print $1}'| head -1")
	   type.command <- paste ("head -200 reference.fasta | ", infoseq.exe)
	   seq.type <- system(type.command, intern = TRUE )
	   if ( seq.type != "P" ){
	   	   stop("CHIPSTER-NOTE: Your reference file is not a fasta fromatted PROTEIN sequence file")
	   }
	   
	   
	   index.command <- paste(diamond.binary, "makedb --in reference.fasta -d reference >> diamond.log")
	   system(index.command)
	   keep_index <- ("yes")
	   
   }	
}

if (use_remote_index == "yes"){

  # make Diamond directry if needed
  if(!dir.exists("/tmp/diamond")){
  	system("mkdir /tmp/diamond")
  }
  if ( db == "swiss" ){
	  dbfile <- ("swiss.dmnd")
  }
  
  if (db == "trembl" ){
	  dbfile <- ("trembl.dmnd")
  }
  
  if ( db == "nr" ){
	  dbfile <- ("nr.dmnd")
  }
  db_in_chipster <- paste( "/tmp/diamond", dbfile, sep="/" )

  #Check if download is needed
  dl.path <- paste(  index_url, dbfile, sep="/" )	
  md5file <- paste(dbfile, "md5", sep=".")
  md5.path <- paste( dl.path, "md5", sep="." )
  if(!file.exists(db_in_chipster)){
  #Download the indexex
	wget.command <- paste("p=$(pwd); cd /tmp/diamond; wget", dl.path, "&>>${p}/diamond.log" )
	system(wget.command)
	wget.command <- paste("p=$(pwd); cd /tmp/diamond; wget", md5.path, "&>>${p}/diamond.log" )
	system(wget.command)
  }else{
	#check md5 file
    wget.command <- paste("wget -O md5test ", md5.path, "&>>${p}/diamond.log" )	
	system(wget.command)
	diff.command <- paste("diff md5test /tmp/diamond/", md5file, " | wc -l ", sep="")
	md5comparison <- system(diff.command,  intern = TRUE )
	if( md5comparison == 0){
		system("echo Database already downloaded >> diamond.log")
		system("ls -l /tmp/diamond >> diamond.log")
		system("rm -f md5test")
	}else{
		system("rm -f md5test")
		rm.command <- paste( "rm -f /tmp/diamond/", dbfile, "*", sep="")
		system(rm.command)
	    wget.command <- paste("p=$(pwd); cd /tmp/diamond; wget", dl.path, "&>>${p}/diamond.log" )	
	    system("echo New version of the database downloaded >> diamond.log")
	    system("ls -l /tmp/diamond >> diamond.log")
    }
}
  
}  

#count the query sequeces
seqcount.exe <- file.path(emboss.path, "seqcount -filter sequence")
str.queryseq <- system(seqcount.exe, intern = TRUE )
num.queryseq <- as.integer(str.queryseq)
#round(num.queryseq)

if (num.queryseq > 1000000){
	stop(paste('CHIPSTER-NOTE: Too many query sequences. Maximun is 1000000 but your file contains ', num.queryseq ))
}


if ( str.filetype != "fasta"){
	seqret.binary <- file.path(emboss.path,"seqret")
	seqret.command <- paste(seqret.binary, "sequence sequence.fasta -auto")
	system(seqret.command)
	system("rm -f sequence")
	system("mv sequence.fasta sequence")
}

#Check if the query is protein or nucleotode
infoseq.exe <- file.path(emboss.path, "infoseq -nowarning -noerror -nodie -only -type -nohead -filter | uniq | awk '{print $1}'| head -1")
type.command <- paste ("head -200 sequence | ", infoseq.exe)
seq.type <- system(type.command, intern = TRUE )

method <- ("blastx")
if ( seq.type == "P" ){
	method <- ("blastp")
}

outfile <- ("diamond.txt")

if ( outformat == "0" ){
	outfile <- ("diamond.txt")
}
if ( outformat == "5" ){
	outfile <- ("diamond.xml")
}
if ( outformat == "6" ){
	outfile <- ("diamond.csv")
}

if ( outformat == "100" ){
	outfile <- ("diamond.txt")
}
if ( outformat == "101" ){
	outfile <- ("diamond.sam")
}

diamond.options <- paste( method, " --query sequence --out ", outfile," --outfmt ", outformat," -p 4 -d ", db, sep="")
if( use_remote_index == "yes"){
	diamond.options <- paste(diamond.options, " -d /tmp/diamond/", db, sep="")
} else {
	diamond.options <- paste(diamond.options, " -d reference")
}
diamond.options <- paste(diamond.options, "--max-target-seqs", maxtargetseqs)
diamond.options <- paste(diamond.options, "--evalue", evalue)
diamond.options <- paste(diamond.options, "--matrix", matrix)
if ( unal == "yes" ){
	diamond.options <- paste(diamond.options, "--unal 1 --un unaligned.fasta ")
}
if (smode == "sensitive"){
	diamond.options <- paste(diamond.options, "--sensitive")
}
if (smode == "moresensitive"){
	diamond.options <- paste(diamond.options, "--more-sensitive")
}




command.full <- paste(diamond.binary, diamond.options, " 2>>diamond.log" )
echo.command <- paste("echo ", diamond.binary, diamond.options, " >>diamond.log" )
system(echo.command)
#command.full <- paste(conda.path, conda.def, " -in sequence ", muscle.options, " 2>>alignment.log" )
system(command.full)
if ( save_log == "no") {
	system ("rm -f diamond.log")
}


if ( keep_index == "no") {
	system ("rm -f reference.dmnd")
}
