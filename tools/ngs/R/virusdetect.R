# TOOL virusdetect.R: "VirusDetect" (VirusDetect pipeline performs virus identification using sRNA sequencing data. Given a FASTQ file, it performs de novo assembly and reference-guided assembly by aligning sRNA reads to the known virus reference database. The assembled contigs are compared to the virus reference sequences with BLAST for virus identification.)
# INPUT inputseq: "Input reads file" TYPE GENERIC (Reads file)
# OUTPUT OPTIONAL virusdetect_contigs.fa 
# OUTPUT OPTIONAL contigs_with_blastn_matches.fa 
# OUTPUT OPTIONAL contigs_with_blastx_matches.fa 
# OUTPUT OPTIONAL undetermined_contigs.fa 
# OUTPUT OPTIONAL blastn_matching_references.fa 
# OUTPUT OPTIONAL blastn_matching_references.fa.fai
# OUTPUT OPTIONAL blastn_matching_references.html 
# OUTPUT OPTIONAL blastx_matching_references.fa 
# OUTPUT OPTIONAL blastx_matching_references.fa.fai
# OUTPUT OPTIONAL blastx_matching_references.html 
# OUTPUT OPTIONAL blastn_matches.bam
# OUTPUT OPTIONAL blastn_matches.bam.bai 
# OUTPUT OPTIONAL blastx_matches.bam 
# OUTPUT OPTIONAL blastx_matches.bam.bai
# OUTPUT OPTIONAL blastn_matches.tsv 
# OUTPUT OPTIONAL blastx_matches.tsv 
# OUTPUT OPTIONAL undetermined.html
# OUTPUT OPTIONAL undetermined_blast.html
# OUTPUT OPTIONAL {...}.pdf
# OUTPUT OPTIONAL vd.log
# OUTPUT OPTIONAL virusdetect_results.tar
# PARAMETER OPTIONAL reference: "Reference virus database" TYPE [vrl_plant: "Plant viruses", vrl_algae: "Algae viruses", vrl_bacteria: "Bacterial viruses", vrl_fungus: "Fungal viruses", vrl_invertebrate: "Invertebrate viruses", vrl_protozoa: "Protozoa viruses", vrl_vertebrate: "Vertebrate viruses"] DEFAULT vrl_plant (Reference virus database.)
# PARAMETER OPTIONAL hostorg: "Host organism" TYPE [none, "FILES genomes/indexes/bwa .fa"] DEFAULT none (Host organism.)
# PARAMETER OPTIONAL hsp_cover: "Minimum fraction of a contig covered by virus reference" TYPE DECIMAL DEFAULT 0.75 (At least this fraction of a contig has to match to the virus reference sequence by BLAST, otherwise the hit is not considered for virus assignment.)
# PARAMETER OPTIONAL coverage_cutoff: "Minimum fraction of virus reference covered by contigs" TYPE DECIMAL DEFAULT 0.1 (At least this fraction of a virus reference sequence has to be covered by contigs, otherwise the virus assignment is not made.)
# PARAMETER OPTIONAL depth_cutoff: "Minimum read depth" TYPE INTEGER DEFAULT 5 (Depth cutoff for virus reference assignment.)  
# PARAMETER OPTIONAL blast_ref: "Return matching reference sequences" TYPE [yes: Yes, no: No] DEFAULT no (Return the reference sequences for BLASTX and BLASTN matches.)
# PARAMETER OPTIONAL blast_bam: "Return BAM formatted alignments" TYPE [yes: Yes, no: No] DEFAULT no (Return the BAM formatted alignments of the contigs to the reference sequences.)
# PARAMETER OPTIONAL save_log: "Collect a log file" TYPE [yes: Yes, no: No] DEFAULT yes (Collect a log file about the analysis run.)
# PARAMETER OPTIONAL sn_tag: "Use input names in output file names" TYPE [yes: Yes, no: No] DEFAULT yes (Name the output files according to the input sequence file.)
# PARAMETER OPTIONAL save_tar: "Return results in one archive file" TYPE [yes: Yes, no: No] DEFAULT yes (Collect all the output files into a single tar formatted file.)


# 13.11.2016 KM, created
# 17.8.2017 KM, option to use input file names for output files

options(scipen=999)

source(file.path(chipster.common.path, "zip-utils.R"))
unzipIfGZipFile("inputseq")

# Handle output names
source(file.path(chipster.common.path, "tool-utils.R"))
# read input names
inputnames <- read_input_definitions()


samtools.binary <- c(file.path(chipster.tools.path, "samtools-0.1.19", "samtools"))
vd.binary <- c(file.path(chipster.tools.path, "virusdetect", "virus_detect.pl"))
#vd.binary <- c(file.path("/opt/chipster/tools_local/virusdetect", "virus_detect.pl"))
vd.parameters <- paste("--reference", reference, "--thread_num", chipster.threads.max, "--coverage_cutoff", coverage_cutoff, "--depth_cutoff", depth_cutoff, "--hsp_cover", hsp_cover )

system("date > vd.log")


#if (save_tar == "yes" ){
#	if (sn_tag == "yes" ){
#	stop("CHIPSTER-NOTE: You can't use tar archive outout formta together with input file name based result file names")
#  }
#}



if (hostorg != "none" ){
	#If host sequence subtraction is used, then we need to create a temporary copy of virus detect
	#vdpath <-  c(file.path("/opt/chipster/tools_local", "virusdetect"))
	#vdpath <-  c(file.path(chipster.tools.path, "virusdetect"))
	#cp.command <- paste("cp -r ",  vdpath , "./")
	#system(cp.command)
	#vd.binary <- c(file.path("./virusdetect", "virus_detect.pl"))
	
	
	# Using pre-calculated indexes
    bwa.genome <- file.path(chipster.tools.path, "genomes", "indexes", "bwa", hostorg)
	echo.command <- paste("echo Using genome:", bwa.genome, " as the host genome. >> vd.log 2>&1" )
	system(echo.command)
		
	bwa.genome.all <- paste(bwa.genome, "*", sep = "" )
	ln.command <- paste("ln -s ", bwa.genome.all, "./ >> vd.log 2>&1 ")
	system(ln.command)
	##hostorg.fa <-  paste("./virusdetect/databases/", hostorg, ".fa", sep = "")
	hostorg.fa <-  paste( hostorg, ".fa", sep = "")
	##mv.command <- paste("mv ", hostorg.fa ," ./virusdetect/databases/", hostorg, sep = "")
	mv.command <- paste("mv ", hostorg.fa ," ", hostorg, sep = "")
	system(mv.command)
	#system("echo Precalculated indexes linked to working directory >> vd.log 2>&1 ")
	##system("ls -l ./virusdetect/databases/ >> vd.log 2>&1 ")
	#system("ls -l >> vd.log 2>&1 ")
	#system("date >> vd.log")
	vd.parameters <- paste(vd.parameters, "--host_reference", hostorg)	
}

vd.parameters <- paste(vd.parameters, "inputseq")

command.full <- paste(vd.binary, vd.parameters, ' >> vd.log 2>&1' )
system("echo Executing VirusDetect command: >> vd.log")
#system("date >> vd.log")
echo.command <- paste('echo "',command.full, ' ">> vd.log' )
system(echo.command)
system(command.full)
system("echo  >> vd.log")
system("echo VirusDetect analysis finished >> vd.log")
system("date >> vd.log")
system("echo  >> vd.log")
system("echo  ---------------------------------------------------- >> vd.log")
system("echo  Following output files were collected:>> vd.log")

nprefix <- ("")
if ( sn_tag == "yes") {
	seq_ifn <- strip_name(inputnames$inputseq)
	nprefix <- paste(seq_ifn, "_", sep= "" )
}

#virusdetect_contigs.fa 
if (file.exists("result_inputseq/contig_sequences.fa")){
	system("mv result_inputseq/contig_sequences.fa  ./virusdetect_contigs.fa ")
	echo.command <- paste("echo ", nprefix, "virusdetect_contigs.fa'\t'Sequences of non-redundant contigs derived through reference-guided and de novo assemblies.>> vd.log 2>&1", sep = "")
	system(echo.command)
}

#contigs_with_blastn_matches.fa
if (file.exists("result_inputseq/contig_sequences.blastn.fa")){
	system("mv result_inputseq/contig_sequences.blastn.fa  ./contigs_with_blastn_matches.fa")
	echo.command <- paste("echo ", nprefix, "contigs_with_blastn_matches.fa'\t'Sequences of contigs that match to virus references by BLASTN. >> vd.log 2>&1", sep = "" )
	system(echo.command)
	
}

#contigs_with_blastx_matches.fa
if (file.exists("result_inputseq/contig_sequences.blastx.fa")){
	system("mv result_inputseq/contig_sequences.blastx.fa  ./contigs_with_blastx_matches.fa")
	echo.command <- paste("echo ", nprefix, "contigs_with_blastx_matches.fa'\t'Sequences of contigs that match to virus references by BLASTX. >> vd.log 2>&1", sep = "" )
	system(echo.command)
}

#undetermined_contigs.fa
if (file.exists("result_inputseq/contig_sequences.undetermined.fa")){
	system("mv result_inputseq/contig_sequences.undetermined.fa  ./undetermined_contigs.fa")
	echo.command <- paste("echo ", nprefix, "undetermined_contigs.fa'\t'Sequences of contigs that do not match to virus references. >> vd.log 2>&1", sep = "" )
	system(echo.command)
}

#blastn_matching_references.html
if (file.exists("result_inputseq/blastn.html")){
	system("echo '<html>' > blastn_matching_references.html")
	system("awk '{ if ( NR > 6 ) print $0 }' result_inputseq/blastn.html >> blastn_matching_references.html")
	system("echo '</html>' >> blastn_matching_references.html")
	echo.command <- paste("echo ", nprefix, "blastn_matching_references.html '\t'Table listing reference viruses that have corresponding virus contigs identified by BLASTN.  >> vd.log 2>&1", sep = "" )
	system(echo.command)
	system("echo '\t\t\t\t'In addition, a pdf formatted report file is returned for each match. >> vd.log 2>&1")
	for.command <- paste("for file in $(ls result_inputseq/blastn_references/*.html); do bln=$(basename $file .html); ", chipster.tools.path, "/miniconda3/envs/chipster_tools/bin/weasyprint $file ${bln}.bn.pdf >> vd.log 2>&1; done;", sep = "")
    system(for.command)
	#system(" for file in $(ls result_inputseq/blastn_references/*.html); do bln=$(basename $file .html); /opt/chipster/tools/miniconda3/envs/chipster_tools/bin/weasyprint $file ${bln}.bn.pdf >> vd.log 2>&1; done;")
}

#blastx_matching_references.html
if (file.exists("result_inputseq/blastx.html")){
	system("echo '<html>' > blastx_matching_references.html")
	system ("awk '{ if ( NR > 6 ) print $0 }' result_inputseq/blastx.html >> blastx_matching_references.html")
	system("echo '</html>' >> blastx_matching_references.html")
	echo.command <- paste("echo ", nprefix, "blastx_matching_references.html '\t'Table listing reference viruses that have corresponding virus contigs identified by BLASTX. >> vd.log 2>&1", sep = "" )
	system(echo.command)
	system("echo '\t\t\t\t'In addition, a pdf formatted report file is returned for each match. >> vd.log 2>&1")
	for.command <- paste(" for file in $(ls result_inputseq/blastx_references/*.html); do bln=$(basename $file .html); ", chipster.tools.path, "/miniconda3/envs/chipster_tools/bin/weasyprint $file ${bln}.bx.pdf; done;", sep="")
	system(for.command)
	#system(" for file in $(ls result_inputseq/blastx_references/*.html); do bln=$(basename $file .html); /opt/chipster/tools/miniconda3/envs/chipster_tools/bin/weasyprint $file ${bln}.bx.pdf; done;")
	
}


if ( blast_ref == "yes") {
	#blastn_matching_references.fa
	if (file.exists("result_inputseq/blastn.reference.fa")){
		system("mv result_inputseq/blastn.reference.fa blastn_matching_references.fa")
		system(paste(samtools.binary, "faidx blastn_matching_references.fa"))
		echo.command <- paste("echo ", nprefix, "blastn_matching_references.fa and .fai.'\t'Virus reference sequences that produced hits for BLASTN search with the potential virus contigs. >> vd.log 2>&1", sep = "" )
		system(echo.command)
	}
	
	
	#blastx.reference.fa
	if (file.exists("result_inputseq/blastx.reference.fa")){
		system("mv result_inputseq/blastx.reference.fa blastx_matching_references.fa")
		system(paste(samtools.binary, "faidx blastx_matching_references.fa"))
		echo.command <- paste("echo ", nprefix, "blastx_matching_references.fa and .fai.'\t'Virus reference sequences that produced hits for BLASTX search with the potential virus contigs. >> vd.log 2>&1", sep = "" )
		system(echo.command)
	}
	
}

if ( blast_bam == "yes") {
	
	samtools.binary <- c(file.path(chipster.tools.path, "samtools-0.1.19", "samtools"))
	#blastn_matches.sam
	if (file.exists("result_inputseq/inputseq.blastn.sam")){
		system ("mv result_inputseq/inputseq.blastn.sam blastn.sam")
		# convert sam to bam
		system(paste(samtools.binary, "view -bS blastn.sam -o blastn.bam"))
		system(paste(samtools.binary, "sort blastn.bam blastn_matches"))
		system(paste(samtools.binary, "index blastn_matches.bam"))
		echo.command <- paste("echo ", nprefix, "blastn_matches.bam and .bai. '\t'BAM file containing the BLASTN alignment of each contig to its corresponding virus reference sequence. >> vd.log 2>&1", sep = "" )
		system(echo.command)	
	}
	
	#blastx_matches.sam
	if (file.exists("result_inputseq/inputseq.blastx.sam")){
		system ("mv result_inputseq/inputseq.blastx.sam blastx.sam")
		# convert sam to bam
		system(paste(samtools.binary, "view -bS blastx.sam -o blastx.bam"))
		system(paste(samtools.binary, "sort blastx.bam blastx_matches"))
		system(paste(samtools.binary, "index blastx_matches.bam"))
		echo.command <- paste("echo ", nprefix, "blastnx_matches.bam and .bai. '\t'BAM file containing the BLASTX alignment of each contig to its corresponding virus reference sequence. >> vd.log 2>&1", sep = "" )
		system(echo.command)
	}
}

#blastn_matches.tsv
if (file.exists("result_inputseq/inputseq.blastn.xls")){
	system("mv result_inputseq/inputseq.blastn.xls blastn_matches.tsv")
	echo.command <- paste("echo ", nprefix, "blastn_matches.tsv '\t\t'Table of BLASTN matches to the reference virus database. >> vd.log 2>&1", sep = "" )
	system(echo.command)
}

#blastx_matches.tsv
if (file.exists("result_inputseq/inputseq.blastx.xls")){
	system("mv result_inputseq/inputseq.blastx.xls blastx_matches.tsv")
	echo.command <- paste("echo ", nprefix, "blastx_matches.tsv '\t\t'Table of BLASTX matches to the reference virus database. >> vd.log 2>&1", sep = "" )
	system(echo.command)
}

#Undetermined
if (file.exists("result_inputseq/undetermined.html")){
	system("echo '<html>' > undetermined.html")
    system("awk '{ if ( NR > 6 ) print $0 }' result_inputseq/undetermined.html >> undetermined.html")
	system("echo '</html>' >> undetermined.html")
	echo.command <- paste("echo ", nprefix, "undetermined.html '\t\t'Table listing the length, siRNA size distribution and 21-22nt percentage of undetermined contigs. >> vd.log 2>&1", sep = "" )
    system(echo.command)
	system("echo  '\t\t\t\t'Potential virus contigs are indicated in green. >> vd.log 2>&1")
}

if (file.exists("result_inputseq/undetermined_blast.html")){	
	system("echo '<html>' > undetermined_blast.html")
	system("awk '{ if ( NR > 6 ) print $0 }' result_inputseq/undetermined_blast.html >> undetermined_blast.html")
	system("echo '</html>' >> undetermined_blast.html")
	echo.command <- paste("echo ", nprefix, "undetermined_blast.html '\t'Table listing contigs having hits in the virus reference database but not assigned to >> vd.log 2>&1", sep = "" )
	system(echo.command)
	system("echo '\t\t\t\t'any reference viruses because they did not meet the coverage or depth criteria. >> vd.log 2>&1")
}


#system ("ls -l >> vd.log")
#system ("ls -l result_inputseq >> vd.log")


if ( save_tar == "yes") {
	if ( sn_tag == "yes") {
		seq_ifn <- strip_name(inputnames$inputseq)
	}
	#system ("echo Collecting results >> vd.log") 
	system ("mkdir vd_output")
	system ("mv virusdetect_contigs.fa vd_output/")
	system ("mv contigs_with_blastn_matches.fa vd_output/")
	system ("mv contigs_with_blastx_matches.fa vd_output/")
	system ("mv undetermined_contigs.fa vd_output/")
	system ("mv blastn_matching_references.fa vd_output/")
	system ("mv blastn_matching_references.fa.fai vd_output/")
	system ("mv blastn_matching_references.html vd_output/")
	system ("mv blastx_matching_references.fa vd_output/")
	system ("mv blastx_matching_references.fa.fai vd_output/")
	system ("mv blastx_matching_references.html vd_output/")
	system ("mv blastn_matches.bam vd_output/")
	system ("mv blastn_matches.bam.bai vd_output/")
	system ("mv blastx_matches.bam vd_output/")
	system ("mv blastx_matches.bam.bai vd_output/")
	system ("mv blastn_matches.tsv vd_output/")
	system ("mv blastx_matches.tsv vd_output/")
	system ("mv undetermined.html vd_output/")
	system ("mv undetermined_blast.html vd_output/")
	system ("mv *.pdf vd_output/")
	#system ("mv vd.log vd_output/")
	seq_ifn <- strip_name(inputnames$inputseq)
	if ( sn_tag == "yes") {
		pdf_name_command <- paste('cd vd_output; for n in *; do mv $n ', seq_ifn, '_$n ; done', sep = "")
		#echo.command <- paste("echo '",pdf_name_command, " '>> vd.log" )
		#system(echo.command)
		system(pdf_name_command)	 
	}
	system ("echo ------------------------------------------------------------ >> vd.log")
	system ("cd vd_output; ls -lh >> ../vd.log")
	system ("cd vd_output; tar cf virusdetect_results.tar ./*; mv virusdetect_results.tar ../virusdetect_results.tar ")
	seq_ifn <- strip_name(inputnames$inputseq)
	# Make a matrix of output names
	outputnames <- matrix(NA, nrow=1, ncol=2)
	outputnames[1,] <- c("virusdetect_results.tar", paste(seq_ifn, "_VD_results.tar", sep =""))
	# Write output definitions file
	write_output_definitions(outputnames)
	system ("echo ------------------------------------------------------------- >> vd.log" )
	system ("echo Results have been collected to a single tar formatted archive file. >> vd.log")
	system ("echo You can use tool: Extract .tar or .tar.gz file in Utilities folder to extract result files from the tar archive. >> vd.log")
	#system ("ls -l >> vd.log")	
	
}

if ( sn_tag == "yes") {
	seq_ifn <- strip_name(inputnames$inputseq)
	#system('ls -l >> vd.log')
	pdf_name_command <- paste('for n in *.pdf; do mv $n ', seq_ifn, '_$n ; done', sep = "")
	#echo.command <- paste("echo '",pdf_name_command, " '>> vd.log" )
	#system(echo.command)
	system(pdf_name_command)
	#system('ls -l >> vd.log')
    # Make a matrix of output names
	outputnames <- matrix(NA, nrow=20, ncol=2)
	outputnames[1,] <- c("virusdetect_contigs.fa", paste(seq_ifn, "virusdetect_contigs.fa", sep ="_"))
	outputnames[2,] <- c("contigs_with_blastn_matches.fa", paste(seq_ifn, "contigs_with_blastn_matches.fa", sep ="_"))
	outputnames[3,] <- c("contigs_with_blastx_matches.fa", paste(seq_ifn, "contigs_with_blastx_matches.fa", sep ="_"))
	outputnames[4,] <- c("undetermined_contigs.fa", paste(seq_ifn, "undetermined_contigs.fa", sep ="_"))
	outputnames[5,] <- c("blastn_matching_references.fa", paste(seq_ifn, "blastn_matching_references.fa", sep ="_"))
	outputnames[6,] <- c("blastn_matching_references.fa.fai", paste(seq_ifn, "blastn_matching_references.fa.fai", sep ="_"))
	outputnames[7,] <- c("blastn_matching_references.html", paste(seq_ifn, "blastn_matching_references.html", sep ="_"))
	outputnames[8,] <- c("blastx_matching_references.fa", paste(seq_ifn, "blastx_matching_references.fa", sep ="_"))
	outputnames[9,] <- c("blastx_matching_references.fa.fai", paste(seq_ifn, "blastx_matching_references.fa.fai", sep ="_"))
	outputnames[10,] <- c("blastx_matching_references.html", paste(seq_ifn, "blastx_matching_references.html", sep ="_"))
	outputnames[11,] <- c("blastn_matches.bam.bai", paste(seq_ifn, "blastn_matches.bam.bai", sep ="_"))
	outputnames[12,] <- c("blastn_matches.bam", paste(seq_ifn, "blastn_matches.bam", sep ="_"))
	outputnames[13,] <- c("blastx_matches.bam", paste(seq_ifn, "blastx_matches.bam", sep ="_"))
	outputnames[14,] <- c("blastx_matches.bam.bai", paste(seq_ifn, "blastx_matches.bam.bai", sep ="_"))
	outputnames[15,] <- c("blastn_matches.tsv", paste(seq_ifn, "blastn_matches.tsv", sep ="_"))
	outputnames[16,] <- c("blastx_matches.tsv", paste(seq_ifn, "blastx_matches.tsv", sep ="_"))
	outputnames[17,] <- c("undetermined.html", paste(seq_ifn, "undetermined.html", sep ="_"))
	outputnames[18,] <- c("undetermined_blast.html", paste(seq_ifn, "undetermined_blast.html", sep ="_"))
    outputnames[19,] <- c("vd.log", paste(seq_ifn, "vd.log", sep ="_"))
    outputnames[20,] <- c("virusdetect_results.tar", paste(seq_ifn, "virusdetect_results.tar", sep ="_"))
	# Write output definitions file
	write_output_definitions(outputnames)	
}

if ( save_log == "no") {
	system ("rm -f vd.log")
}

