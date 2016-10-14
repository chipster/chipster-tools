# TOOL virusdetect.R: "VirusDetect" (VirusDetect  analyzes large-scale sRNA datasets for virus identification. The program performs reference-guided assembly by aligning sRNA reads to the known virus reference database as well as de novo assembly. The assembled contigs are compared to the reference virus sequences for virus identification.  )
# INPUT inputseq: "Input reads file" TYPE GENERIC (Reads file)
# INPUT OPTIONAL hostgenome: "Host genome, optional " TYPE GENERIC (Host genome used for host sequence subtraction )
# OUTPUT OPTIONAL virusdetect_contigs.fa 
# OUTPUT OPTIONAL virusderect_matches_blastn.fa 
# OUTPUT OPTIONAL virusderect_matches_blastx.fa 
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
# PARAMETER OPTIONAL reference: "Reference virus database" TYPE [vrl_plant: "Plant viruses", vrl_algae: "Algae viruses", vrl_bacteria: "Bacterial viruses", vrl_fungus: "Fungal viruses", vrl_invertebrate: "Invertebrate viruses", vrl_protozoa: "Protozoa viruses", vrl_vertebrate: "Vertebrate viruses"] DEFAULT vrl_plant (Reference virus database.)
# PARAMETER OPTIONAL hostorg:  "Host organism" TYPE [none, other, Arabidopsis_thaliana.TAIR10.30, Bos_taurus.UMD3.1, Canis_familiaris.BROADD2.67, Canis_familiaris.CanFam3.1, Drosophila_melanogaster.BDGP5, Drosophila_melanogaster.BDGP6, Felis_catus.Felis_catus_6.2, Gallus_gallus.Galgal4, Gasterosteus_aculeatus.BROADS1, Halorubrum_lacusprofundi_atcc_49239.GCA_000022205.1.30, Homo_sapiens.GRCh37.75, Homo_sapiens.GRCh38, Homo_sapiens.NCBI36.54, mature, Medicago_truncatula.GCA_000219495.2.30, Mus_musculus.GRCm38, Mus_musculus.NCBIM37.67, Oryza_sativa.IRGSP-1.0.30, Ovis_aries.Oar_v3.1, Populus_trichocarpa.JGI2.0.30, Rattus_norvegicus.RGSC3.4.69, Rattus_norvegicus.Rnor_5.0, Rattus_norvegicus.Rnor_6.0, Schizosaccharomyces_pombe.ASM294v2.30, Solanum_tuberosum.3.0.30, Sus_scrofa.Sscrofa10.2, Vitis_vinifera.IGGP_12x.30, Yersinia_enterocolitica_subsp_palearctica_y11.GCA_000253175.1.30, Yersinia_pseudotuberculosis_ip_32953_gca_000834295.GCA_000834295.1.30] DEFAULT none (Reference sequence.)
# PARAMETER OPTIONAL hsp_cover: "Reference virus coverage cuttoff" TYPE DECIMAL DEFAULT 0.75 (Coverage cutoff of a reported virus contig by reference virus sequences.)
# PARAMETER OPTIONAL coverage_cutoff: "Assembled virus contig cuttoff" TYPE DECIMAL DEFAULT 0.1 (Coverage cutoff of a reported virus reference sequences by assembled virus contigs. )
# PARAMETER OPTIONAL depth_cutoff: "Depth cutoff" TYPE INTEGER DEFAULT 5 (Depth cutoff of a reported virus reference)  
# PARAMETER OPTIONAL blast_ref: "return the matching reference sequences" TYPE [yes: Yes, no: No] DEFAULT no (Return the reference sequences for BLASTx and BLASTn runs)
# PARAMETER OPTIONAL blast_bam: "return the BAM formattet alignments" TYPE [yes: Yes, no: No] DEFAULT no (Return the BAM formatted alignments of the viral sequeces to the refernce sequences)
# PARAMETER OPTIONAL save_log: "Collect a log file" TYPE [yes: Yes, no: No] DEFAULT no (Collect a log file about the analysis run.)

# KM 2.9. 2016

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
vd.parameters <- paste("--reference", reference, "--thread-num", chipster.threads.max, "--coverage-cutoff", coverage_cutoff, "--depth-cutoff", depth_cutoff, "--hsp-cover", hsp_cover )

system("date > vd.log")


if (hostorg == "none" ){
	if (file.exists("hostgenome")){
		stop(paste('CHIPSTER-NOTE: ', "You have defined a host genome file but the host organism genome setting has value: none."))
	}
	
}


if (hostorg != "none" ){
	#If host sequence subtraction is used, then we need to cerate a temporary copy of virus detect
	#vdpath <-  c(file.path("/opt/chipster/tools_local", "virusdetect"))
	#vdpath <-  c(file.path(chipster.tools.path, "virusdetect"))
	#cp.command <- paste("cp -r ",  vdpath , "./")
	#system(cp.command)
	#vd.binary <- c(file.path("./virusdetect", "virus_detect.pl"))
	
	if (hostorg == "other"){
		if (file.exists("hostgenome")){
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
		}else{
			stop(paste('CHIPSTER-NOTE: ', "You need to provide a FASTA file or choose one of the provided reference genomes."))
		} 
		
	}else{
		# Using pre-calculated indexes
		bwa.genome <- file.path(chipster.tools.path, "genomes", "indexes", "bwa", hostorg)
		echo.command <- paste("echo ", bwa.genome, " >> vd.log 2>&1" )
		system(echo.command)
		
		bwa.genome.all <- paste(bwa.genome, "*", sep = "" )
		#ls.command <- paste("ls -l ", bwa.genome.all, ">> vd.log 2>&1 " )
		#system(ls.command)
		
		#cp.command <- paste("cp", bwa.genome.all, "./virusdetect/databases/ >> vd.log 2>&1 ")
		cp.command <- paste("cp", bwa.genome.all, "./ >> vd.log 2>&1 ")
		system(cp.command)
		#hostorg.fa <-  paste("./virusdetect/databases/", hostorg, ".fa", sep = "")
		hostorg.fa <-  paste( hostorg, ".fa", sep = "")
		#mv.command <- paste("mv ", hostorg.fa ," ./virusdetect/databases/", hostorg, sep = "")
		mv.command <- paste("mv ", hostorg.fa ," ", hostorg, sep = "")
		system(mv.command)
		system("echo Precalculated indexes copied to working directory >> vd.log 2>&1 ")
		#system("ls -l ./virusdetect/databases/ >> vd.log 2>&1 ")
		#system("ls -l >> vd.log 2>&1 ")
		system("date >> vd.log")
		vd.parameters <- paste(vd.parameters, "--host-reference", hostorg)
	}
	
	
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





