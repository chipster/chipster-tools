# TOOL virusdetect-with-index-building.R: "VirusDetect with own host genome" (VirusDetect pipeline performs virus identification using sRNA sequencing data. Given a FASTQ file, it performs de novo assembly and reference-guided assembly by aligning sRNA reads to the known virus reference database. The assembled contigs are compared to the reference virus sequences for virus identification.)
# INPUT inputseq: "Input reads file" TYPE GENERIC (Reads file)
# INPUT hostgenome: "Host genome " TYPE GENERIC (Host genome used for host sequence subtraction. This can be a fasta formatted sequence file or a BWA index file created by Chipster.)
# OUTPUT OPTIONAL virusdetect_contigs.fa 
# OUTPUT OPTIONAL virusdetect_matches_blastn.fa 
# OUTPUT OPTIONAL virusdetect_matches_blastx.fa 
# OUTPUT OPTIONAL contig_sequences.undetermined.fa 
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
# OUTPUT OPTIONAL hostgenome_bwa_index.tar
# OUTPUT OPTIONAL {...}.pdf
# OUTPUT OPTIONAL vd.log
# OUTPUT OPTIONAL virusdetect_results.tar
# PARAMETER OPTIONAL reference: "Reference virus database" TYPE [vrl_plant: "Plant viruses", vrl_algae: "Algae viruses", vrl_bacteria: "Bacterial viruses", vrl_fungus: "Fungal viruses", vrl_invertebrate: "Invertebrate viruses", vrl_protozoa: "Protozoa viruses", vrl_vertebrate: "Vertebrate viruses"] DEFAULT vrl_plant (Reference virus database.)
# PARAMETER OPTIONAL hsp_cover: "Reference virus coverage cuttoff" TYPE DECIMAL DEFAULT 0.75 (Coverage cutoff of a reported virus contig by reference virus sequences.)
# PARAMETER OPTIONAL coverage_cutoff: "Assembled virus contig cuttoff" TYPE DECIMAL DEFAULT 0.1 (Coverage cutoff of a reported virus reference sequence by assembled virus contigs.)
# PARAMETER OPTIONAL depth_cutoff: "Depth cutoff" TYPE INTEGER DEFAULT 5 (Depth cutoff of a reported virus reference.)  
# PARAMETER OPTIONAL blast_ref: "Return matching reference sequences" TYPE [yes: Yes, no: No] DEFAULT no (Return the reference sequences for BLASTx and BLASTn runs.)
# PARAMETER OPTIONAL blast_bam: "Return BAM formatted alignments" TYPE [yes: Yes, no: No] DEFAULT no (Return the BAM formatted alignments of the viral sequences to the reference sequences.)
# PARAMETER OPTIONAL save_log: "Collect a log file" TYPE [yes: Yes, no: No] DEFAULT no (Collect a log file about the analysis run.)
# PARAMETER OPTIONAL sn_tag: "Use input names in output file names" TYPE [yes: Yes, no: No] DEFAULT yes (Name the output files according to the input sequence file.)
# PARAMETER OPTIONAL save_tar: "Return results in one archive file" TYPE [yes: Yes, no: No] DEFAULT no (Collect all the output into a single tar formatted file.)

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

#check that output options don't conflict
#if (save_tar == "yes" ){
#	if (sn_tag == "yes" ){
#		stop("CHIPSTER-NOTE: You can't use tar archive outout formta together with input file name based result file names")
#	}
#}

bwa.index.binary <- file.path(chipster.module.path, "shell", "check_bwa_index.sh")
unzipIfGZipFile("hostgenome")
hostgenome.filetype <- system("file -b hostgenome | cut -d ' ' -f2", intern = TRUE )
			
# case 1. Ready calculated indexes in tar format
if (hostgenome.filetype == "tar"){
	check.command <- paste ( bwa.index.binary, "hostgenome| tail -1 ")
	bwa.genome <- system(check.command, intern = TRUE)
				
# case 2. Fasta file
}else{
	check.command <- paste ( bwa.index.binary, "hostgenome -tar| tail -1 ")
	bwa.genome <- system(check.command, intern = TRUE)
	cp.command <- paste("cp ", bwa.genome, "_bwa_index.tar ./hostgenome_bwa_index.tar ", sep ="")
	system(cp.command)
	hg_ifn <- strip_name(inputnames$hostgenome)
	# Make a matrix of output names
	outputnames <- matrix(NA, nrow=1, ncol=2)
	outputnames[1,] <- c("hostgenome_bwa_index.tar", paste(hg_ifn, "_bwa_index.tar", sep =""))
	# Write output definitions file
	write_output_definitions(outputnames)
}
			
vd.parameters <- paste(vd.parameters, "--host-reference hostgenome")
#system("ls -l >> vd.log")
system("date >> vd.log")


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

#virusderect_matches_blastn.fa
if (file.exists("result_inputseq/contig_sequences.blastn.fa")){
	system("mv result_inputseq/contig_sequences.blastn.fa  ./virusderect_matches_blastn.fa")
}

#virusderect_matches_blastx.fa
if (file.exists("result_inputseq/contig_sequences.blastx.fa")){
	system("mv result_inputseq/contig_sequences.blastx.fa  ./virusderect_matches_blastx.fa")
}

#contig_sequences.undetermined.fa
if (file.exists("result_inputseq/contig_sequences.undetermined.fa")){
	system("mv result_inputseq/contig_sequences.undetermined.fa  ./contig_sequences.undetermined.fa")
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
	system("mv result_inputseq/inputseq.blastx.tsv blastx_matches.tsv")
}

#system ("ls -l >> vd.log")
#system ("ls -l result_inputseq >> vd.log")

if ( save_log == "no") {
	system ("rm -f vd.log")
}

if ( save_tar == "yes") {
	system ("echo Collecting results >> vd.log") 
	system ("ls -l >> vd.log")
	system ("mkdir vd_output")
	system ("mv virusdetect_contigs.fa vd_output/")
	system ("mv virusderect_matches_blastn.fa vd_output/")
	system ("mv virusderect_matches_blastx.fa vd_output/")
	system ("mv contig_sequences.undetermined.fa vd_output/")
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
	system ("mv hostgenome_bwa_index.tar vd_output/")
	system ("mv *.pdf vd_output/")
	system ("cp vd.log vd_output/")
	
	seq_ifn <- strip_name(inputnames$inputseq)
	if ( sn_tag == "yes") {
		pdf_name_command <- paste('cd vd_output; for n in *; do mv $n ', seq_ifn, '_$n ; done', sep = "")
		echo.command <- paste("echo '",pdf_name_command, " '>> vd.log" )
		system(echo.command)
		system(pdf_name_command)	 
	}
	system ("cd vd_output; tar cf virusdetect_results.tar ./*; mv virusdetect_results.tar ../virusdetect_results.tar ")
	
	outputnames <- matrix(NA, nrow=1, ncol=2)
	outputnames[1,] <- c("virusdetect_results.tar", paste(seq_ifn, "_VD_results.tar", sep =""))
	# Write output definitions file
	write_output_definitions(outputnames)
	system ("echo Result collectin ready >> vd.log")
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
	outputnames <- matrix(NA, nrow=17, ncol=2)
	outputnames[1,] <- c("virusdetect_contigs.fa", paste(seq_ifn, "virusdetect_contigs.fa", sep ="_"))
	outputnames[2,] <- c("virusderect_matches_blastn.fa", paste(seq_ifn, "virusderect_matches_blastn.fa", sep ="_"))
	outputnames[3,] <- c("virusderect_matches_blastx.fa", paste(seq_ifn, "virusderect_matches_blastx.fa", sep ="_"))
	outputnames[4,] <- c("contig_sequences.undetermined.fa", paste(seq_ifn, "virusderect_matches_blastx.fa", sep ="_"))
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
	outputnames[17,] <- c("vd.log", paste(seq_ifn, "vd.log", sep ="_"))
	# Write output definitions file
	write_output_definitions(outputnames)	
}




