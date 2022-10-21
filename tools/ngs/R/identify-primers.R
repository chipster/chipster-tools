# TOOL identify-primers.R: "Identify primers and the correct orientation" (Given a FASTQ file and the adapter or primer sequences, this tool checks if the adapter/primer sequences were found in the sample and in which orientation. The FASTQ files can be in a tar package but then you need to specify if the reads are paired end or not. )
# INPUT OPTIONAL reads.tar: "Tar package containing the FASTQ files" TYPE GENERIC
# INPUT OPTIONAL file1.fastq: "Forward fastq file" TYPE GENERIC ("Forward fastq file")
# INPUT OPTIONAL file2.fastq: "Reverse fastq file" TYPE GENERIC ("Reverse fastq file")
# OUTPUT primer_summary.tsv
# OUTPUT orientations_summary.txt
# PARAMETER paired: "Is the data paired end or single end reads" TYPE [paired, single] DEFAULT single (Are all the reads paired end? If yes, then identify primers from forward and reverse reads.)
# PARAMETER OPTIONAL adapter5: "The 5' adapter:" TYPE STRING (Give here the 5´ end adapter/primer)
# PARAMETER OPTIONAL adapter3: "The 3' adapter:" TYPE STRING (Give here the 3´ end adapter/primer)
# RUNTIME R-4.1.1-asv

# ES 13.10.2022
# could use either tar package or fastq files as input, Decided to use tar, because it's more convenient in the pipeline
# PARAMETER paired: "Is the data paired end or single end reads" TYPE [paired, single] DEFAULT single (Are all the reads paired end, so one forward and one reverse FASTQ file for one sample. If single end reads,use only those forward parameters.)

source(file.path(chipster.common.path,"tool-utils.R"))
source(file.path(chipster.common.path,"zip-utils.R"))

# load 2 libraries
library(Biostrings)
library(ShortRead)

if (adapter3=="" && adapter5==""){
  stop(paste('CHIPSTER-NOTE: ',"You need to give at least one adapter/primer sequence as a parameter"))
}
#check out if the file is compressed and if so unzip it


reverse_file=""

# if tar package given, use that
if (fileOk("reads.tar")){
    unzipIfGZipFile("reads.tar")
    # Make an input  folder 
    system("mkdir input_folder")

    #untar the tar package to input_folder and list the filenames
    untar("reads.tar", exdir = "input_folder")

    #list the full file names with the relative path / every second file forward, reverse....
    filenames <- sort(list.files("input_folder", full.names=TRUE))
    forward_file = filenames[1]

    if (paired=="paired"){
        reverse_file=filenames[2]
    }    
}else if (fileOk("file1.fastq")){
    unzipIfGZipFile("file1.fastq")
    forward_file <- "file1.fastq"
    if (fileOk("file2.fastq")){
        reverse_file <- "file2.fastq"
    }
}else{
    stop(paste('CHIPSTER-NOTE: ',"No input file recognized. Please give either a tar package or the forward fastq file as input"))
}

# to ensure we have the right primers, and the correct orientation of the primers 
allOrients <- function(primer) {
    # Create all orientations of the input sequence
    require(Biostrings)
    dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
    orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
        RevComp = reverseComplement(dna))
    return(sapply(orients, toString))  # Convert back to character vector
}
# make a summary txt file for all the different orientations. User can copy these and use in next steps
sink("orientations_summary.txt")
    cat("Check all orientations of the given primers/adapters.\n\n")
    if (paired=="paired"){
        cat("Reads are paired end. \n")
        cat("The forward file used to check the adapters:",basename(forward_file),"\n")
        cat("The reverse file used to check the adapters:",basename(reverse_file),"\n\n")
    }else{
        cat("Reads are single end. \n")
        cat("The forward file used to check the adapters:",basename(forward_file),"\n\n")
    }
    if (adapter5!=""){
        cat("All orientations of the 5' adapter:\n")
        orients5 <-allOrients(adapter5)
        cat("Forward: ",orients5[1], "\n")
        cat("Complement: ",orients5[2], "\n")
        cat("Reverse: ",orients5[3], "\n")
        cat("RevComp: ",orients5[4], "\n\n")
        
    }
    if ((adapter3!="")){
        cat("\nAll orientations of the 3' adapter:\n")
        orients3 <-allOrients(adapter3)
        cat("Forward: ",orients3[1], "\n")
        cat("Complement: ",orients3[2], "\n")
        cat("Reverse: ",orients3[3], "\n")
        cat("RevComp: ",orients3[4], "\n\n")
    }
   
sink()

primerHits <- function(primer, fn) {

    # Counts number of reads in which the primer is found
    nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
    return(sum(nhits > 0))
}


# make a tsv table how many hits were found and in which orientation
if (reverse_file==""){

     #single reads so no reverse read here to specify
    if (adapter5==""){
        #Forward_reads_3adapter= t(sapply(orients3, primerHits, fn = forward_file))
        #rownames(Forward_reads_3adapter)<-"3xx"
        #write.table(Forward_reads_3adapter, file ="primer_summary.tsv", sep='\t')
        matrix <- t(data.frame(Forwardreads_3adapter= sapply(orients3, primerHits, fn = forward_file)))
        #names(matrix) <- c("Forward","Complement", "Reverse", "RevComp")
    }else if (adapter3==""){
        matrix <- rbind(Forward_reads_5adapter = sapply(orients5, primerHits, fn = forward_file))
       
    }else{
      
        matrix <- rbind(Forward_reads_5adapter = sapply(orients5, primerHits, fn = forward_file),
                        Forward_reads_3adapter = sapply(orients3, primerHits, fn = forward_file))
    }
}else{  #paired reads, reverse and forward read to check
   
    if (adapter5==""){
        
        matrix <-rbind(Forward_reads.3adapter = sapply(orients3, primerHits,fn= forward_file),
                Reverse_reads_3adapter = sapply(orients3, primerHits, fn = reverse_file))
    
    }else if (adapter3==""){
      
        matrix <-rbind(Forward_reads.5adapter = sapply(orients5, primerHits, fn = forward_file),
            Reverse_reads_5adapter = sapply(orients5, primerHits, fn = reverse_file))
    }else{
        
        matrix <-rbind(Forward_reads_5adapter = sapply(orients5, primerHits, fn = forward_file),
                Forward_reads_3adapter = sapply(orients3, primerHits, fn = forward_file),
                Reverse_reads_5adapter = sapply(orients5, primerHits, fn = reverse_file),
                Reverse_reads_3adapter = sapply(orients3, primerHits, fn = reverse_file))
    }

}

write.table(matrix, file ="primer_summary.tsv", sep='\t')

#EOF