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
# PARAMETER OPTIONAL hsp_cover: "Minimum fraction of a contig covered by reference" TYPE DECIMAL DEFAULT 0.75 (At least this fraction of a contig has to match to the virus reference sequence by BLAST, otherwise the hit is not considered for virus assignment.)
# PARAMETER OPTIONAL coverage_cutoff: "Minimum fraction of reference covered by contigs" TYPE DECIMAL DEFAULT 0.1 (At least this fraction of a virus reference sequence has to be covered by contigs, otherwise the virus assignment is not made.)
# PARAMETER OPTIONAL depth_cutoff: "Depth cutoff" TYPE INTEGER DEFAULT 5 (Normalized depth cutoff for virus reference assignment.)  
# PARAMETER OPTIONAL blast_ref: "Return matching reference sequences" TYPE [yes: Yes, no: No] DEFAULT no (Return the reference sequences for BLASTX and BLASTN matches.)
# PARAMETER OPTIONAL blast_bam: "Return BAM formatted alignments" TYPE [yes: Yes, no: No] DEFAULT no (Return the BAM formatted alignments of the contigs to the reference sequences.)
# PARAMETER OPTIONAL save_log: "Collect a log file" TYPE [yes: Yes, no: No] DEFAULT no (Collect a log file about the analysis run.)
# PARAMETER OPTIONAL sn_tag: "Use input names in output file names" TYPE [yes: Yes, no: No] DEFAULT yes (Name the output files according to the input sequence file.)
# PARAMETER OPTIONAL save_tar: "Return results in one archive file" TYPE [yes: Yes, no: No] DEFAULT no (Collect all the output files into a single tar formatted file.)


# 13.11.2016 KM, created
# 17.8.2017 KM, option to use input file names for output files

options(scipen=999)

source(file.path(chipster.common.path, "zip-utils.R"))
unzipIfGZipFile("inputseq")

# Handle output names
source(file.path(chipster.common.path, "tool-utils.R"))
# read input names
inputnames <- read_input_definitions()


samtools.binary <- c(file.path(chipster.tools.path, "samtools", "samtools"))
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
	echo.command <- paste("echo Using genome:", bwa.genome, " >> vd.log 2>&1" )
	system(echo.command)
		
	bwa.genome.all <- paste(bwa.genome, "*", sep = "" )
	ln.command <- paste("ln -s ", bwa.genome.all, "./ >> vd.log 2>&1 ")
	system(ln.command)
	##hostorg.fa <-  paste("./virusdetect/databases/", hostorg, ".fa", sep = "")
	hostorg.fa <-  paste( hostorg, ".fa", sep = "")
	##mv.command <- paste("mv ", hostorg.fa ," ./virusdetect/databases/", hostorg, sep = "")
	mv.command <- paste("mv ", hostorg.fa ," ", hostorg, sep = "")
	system(mv.command)
	system("echo Precalculated indexes linked to working directory >> vd.log 2>&1 ")
	##system("ls -l ./virusdetect/databases/ >> vd.log 2>&1 ")
	system("ls -l >> vd.log 2>&1 ")
	system("date >> vd.log")
	vd.parameters <- paste(vd.parameters, "--host_reference", hostorg)	
}

vd.parameters <- paste(vd.parameters, "inputseq")

command.full <- paste(vd.binary, vd.parameters, ' >> vd.log 2>&1' )
system("echo Starting virus detect. >> vd.log")
system("date >> vd.log")
echo.command <- paste('echo "',command.full, ' ">> vd.log' )
system(echo.command)
system(command.full)


#virusdetect_contigs.fa 
if (file.exists("result_inputseq/contig_sequences.fa")){
	system("mv result_inputseq/contig_sequences.fa  ./virusdetect_contigs.fa ")
}

#contigs_with_blastn_matches.fa
if (file.exists("result_inputseq/contig_sequences.blastn.fa")){
	system("mv result_inputseq/contig_sequences.blastn.fa  ./contigs_with_blastn_matches.fa")
}

#contigs_with_blastx_matches.fa
if (file.exists("result_inputseq/contig_sequences.blastx.fa")){
	system("mv result_inputseq/contig_sequences.blastx.fa  ./contigs_with_blastx_matches.fa")
}

#undetermined_contigs.fa
if (file.exists("result_inputseq/contig_sequences.undetermined.fa")){
	system("mv result_inputseq/contig_sequences.undetermined.fa  ./undetermined_contigs.fa")
}

#blastn_matching_references.html
if (file.exists("result_inputseq/blastn.html")){
	system("echo '<html>' > blastn_matching_references.html")
	system("awk '{ if ( NR > 6 ) print $0 }' result_inputseq/blastn.html >> blastn_matching_references.html")
	system("echo '</html>' >> blastn_matching_references.html")
	system(" for file in $(ls result_inputseq/blastn_references/*.html); do bln=$(basename $file .html); weasyprint $file ${bln}.bn.pdf; done;")
}

#blastx_matching_references.html
if (file.exists("result_inputseq/blastx.html")){
	system("echo '<html>' > blastx_matching_references.html")
	system ("awk '{ if ( NR > 6 ) print $0 }' result_inputseq/blastx.html >> blastx_matching_references.html")
	system("echo '</html>' >> blastx_matching_references.html")
	system(" for file in $(ls result_inputseq/blastx_references/*.html); do bln=$(basename $file .html); weasyprint $file ${bln}.bx.pdf; done;")
	
}


if ( blast_ref == "yes") {
	#blastn_matching_references.fa
	if (file.exists("result_inputseq/blastn.reference.fa")){
		system("mv result_inputseq/blastn.reference.fa blastn_matching_references.fa")
		system(paste(samtools.binary, "faidx blastn_matching_references.fa"))
	}
	
	
	#blastx.reference.fa
	if (file.exists("result_inputseq/blastx.reference.fa")){
		system("mv result_inputseq/blastx.reference.fa blastx_matching_references.fa")
		system(paste(samtools.binary, "faidx blastx_matching_references.fa"))
	}
	
}

if ( blast_bam == "yes") {
	
	samtools.binary <- c(file.path(chipster.tools.path, "samtools", "samtools"))
	#blastn_matches.sam
	if (file.exists("result_inputseq/inputseq.blastn.sam")){
		system ("mv result_inputseq/inputseq.blastn.sam blastn.sam")
		# convert sam to bam
		system(paste(samtools.binary, "view -bS blastn.sam -o blastn.bam"))
		system(paste(samtools.binary, "sort blastn.bam blastn_matches"))
		system(paste(samtools.binary, "index blastn_matches.bam"))
	}
	
	#blastx_matches.sam
	if (file.exists("result_inputseq/inputseq.blastx.sam")){
		system ("mv result_inputseq/inputseq.blastx.sam blastx.sam")
		# convert sam to bam
		system(paste(samtools.binary, "view -bS blastx.sam -o blastx.bam"))
		system(paste(samtools.binary, "sort blastx.bam blastx_matches"))
		system(paste(samtools.binary, "index blastx_matches.bam"))
	}
}

#blastn_matches.tsv
if (file.exists("result_inputseq/inputseq.blastn.xls")){
	system("mv result_inputseq/inputseq.blastn.xls blastn_matches.tsv")
}

#blastx_matches.tsv
if (file.exists("result_inputseq/inputseq.blastx.xls")){
	system("mv result_inputseq/inputseq.blastx.xls blastx_matches.tsv")
}

#Undetermined
if (file.exists("result_inputseq/undetermined.html")){
	system("echo '<html>' > undetermined.html")
    system("awk '{ if ( NR > 1 ) print $0 }' result_inputseq/undetermined.html >> undetermined.html")
	system("echo '</html>' >> undetermined.html")
}


if (file.exists("result_inputseq/undetermined_blast.html")){
	
	system("echo '<html>' > undetermined_blast.html")
	system("awk '{ if ( NR > 1 ) print $0 }' result_inputseq/undetermined_blast.html >> undetermined_blast.html")
	system("echo '</html>' >> undetermined_blast.html")
}


#system ("ls -l >> vd.log")
#system ("ls -l result_inputseq >> vd.log")


if ( save_tar == "yes") {
	if ( sn_tag == "yes") {
		seq_ifn <- strip_name(inputnames$inputseq)
	}
	system ("echo Collecting results >> vd.log") 
	system ("ls -l >> vd.log")
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
		echo.command <- paste("echo '",pdf_name_command, " '>> vd.log" )
		system(echo.command)
		system(pdf_name_command)	 
	}
	system ("cd vd_output; tar cf virusdetect_results.tar ./*; mv virusdetect_results.tar ../virusdetect_results.tar ")
	seq_ifn <- strip_name(inputnames$inputseq)
	# Make a matrix of output names
	outputnames <- matrix(NA, nrow=1, ncol=2)
	outputnames[1,] <- c("virusdetect_results.tar", paste(seq_ifn, "_VD_results.tar", sep =""))
	# Write output definitions file
	write_output_definitions(outputnames)
	system ("echo Result collecti ready >> vd.log")
	system ("ls -l >> vd.log")	
	
}

if ( sn_tag == "yes") {
	seq_ifn <- strip_name(inputnames$inputseq)
	system('ls -l >> vd.log')
	pdf_name_command <- paste('for n in *.pdf; do mv $n ', seq_ifn, '_$n ; done', sep = "")
	echo.command <- paste("echo '",pdf_name_command, " '>> vd.log" )
	system(echo.command)
	system(pdf_name_command)
	system('ls -l >> vd.log')
    # Make a matrix of output names
	outputnames <- matrix(NA, nrow=19, ncol=2)
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

	if ( save_log == "yes") {
		outputnames[19,] <- c("vd.log", paste(seq_ifn, "vd.log", sep ="_"))
	}
	# Write output definitions file
	write_output_definitions(outputnames)	
}

if ( save_log == "no") {
	system ("rm -f vd.log")
}

