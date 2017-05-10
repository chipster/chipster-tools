# TOOL picard-merge-bam-alignment.R: "Merge BAM alignment" (Merge sorted BAM alignment and unaligned, tagged BAM file. Make sure the input files are assigned correctly!)
# INPUT unmapped.bam: "Unmapped BAM" TYPE GENERIC
# INPUT aligned.bam: "Aligned BAM" TYPE GENERIC
# OUTPUT OPTIONAL merged.bam     
# PARAMETER OPTIONAL reference: "Reference" TYPE [Homo_sapiens.GRCh38, Mus_musculus.GRCm38] DEFAULT Homo_sapiens.GRCh38 (Use same reference as in the alignment!)


# OUTPUT OPTIONAL log.txt

# 2016-10-31 ML
# 2017-05-08 ML add sort BAM to this tool

## NOTE!!! Add the reference genomes!!! Like this:
# PARAMETER organism: "Genome" TYPE [Arabidopsis_thaliana.TAIR10.30, Bos_taurus.UMD3.1.83, Canis_familiaris.CanFam3.1.83, Drosophila_melanogaster.BDGP6.83, Felis_catus.Felis_catus_6.2.83, Gallus_gallus.Galgal4.83, Gasterosteus_aculeatus.BROADS1.83, Halorubrum_lacusprofundi_atcc_49239.GCA_000022205.1.30, Homo_sapiens.GRCh37.75, Homo_sapiens.GRCh38.83, Medicago_truncatula.GCA_000219495.2.30, Mus_musculus.GRCm38.83, Oryza_sativa.IRGSP-1.0.30, Ovis_aries.Oar_v3.1.83, Populus_trichocarpa.JGI2.0.30, Rattus_norvegicus.Rnor_5.0.79, Rattus_norvegicus.Rnor_6.0.83, Schizosaccharomyces_pombe.ASM294v2.30, Solanum_tuberosum.3.0.30, Sus_scrofa.Sscrofa10.2.83, Vitis_vinifera.IGGP_12x.30, Yersinia_enterocolitica_subsp_palearctica_y11.GCA_000253175.1.30, Yersinia_pseudotuberculosis_ip_32953_gca_000834295.GCA_000834295.1.30] DEFAULT Homo_sapiens.GRCh38.83 (Genome or transcriptome that you would like to align your reads against.)


picard.binary <- file.path(chipster.tools.path, "picard-tools", "picard.jar")
path.reference <- "/opt/chipster/genomes/fasta" #/Homo_sapiens.GRCh38.fa"
#path.reference <- "/opt/chipster/genomes/fasta/Homo_sapiens.GRCh38.fa"
#path.tophat.index <- c(file.path(chipster.tools.path, "genomes", "indexes", "tophat2", organism))

# create symlink

#command <- paste("ln -s ", path.reference, "/Homo_sapiens.GRCh38.fa Homo_sapiens.GRCh38.fa", sep="" )
command <- paste("ln -s ", path.reference, "/", reference, ".fa ", reference, ".fa", sep="" )
system(command)		

# create dictionary (dict) (takes less than minute)
#command <- paste("java -Xmx2g -jar ", picard.binary, " CreateSequenceDictionary R=Homo_sapiens.GRCh38.fa O=Homo_sapiens.GRCh38.dict 2>> log.txt", sep="")
command <- paste("java -Xmx2g -jar ", picard.binary, " CreateSequenceDictionary R=", reference, ".fa O=" ,reference, ".dict 2>> log.txt", sep="")
system(command)


# Sort BAM (Picard):
command <- paste("java -Xmx2g -jar", picard.binary, "SortSam I=aligned.bam O=aligned_sorted.bam SO=queryname 2>> log.txt")
system(command)


# Merge files (Picard):
command <- paste("java -Xmx2g -jar ", picard.binary, " MergeBamAlignment UNMAPPED_BAM=unmapped.bam ALIGNED_BAM=aligned_sorted.bam O=merged.bam R=" ,reference, ".fa INCLUDE_SECONDARY_ALIGNMENTS=false PAIRED_RUN=false 2>> log.txt", sep="") 
system(command)
# stop(paste('CHIPSTER-NOTE: ', command))


# EOF