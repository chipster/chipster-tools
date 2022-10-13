# TOOL identify-primers.R: "Identify primers and the correct orientation" (Given a fastg file, this tool check if the given primer was in the .)
# INPUT file1.fastq: "Forward fastq file" TYPE GENERIC ("Forward fastq file")
# INPUT OPTIONAL file2.fastq: "Reverse fastq file" TYPE GENERIC ("Reverse fastq file")
# OUTPUT primer_summary.tsv
# OUTPUT orientations_summary.txt
# PARAMETER OPTIONAL adapter5: "Forward or the 5' adapter to be trimmed" TYPE STRING
# PARAMETER OPTIONAL adapter3: "The 3' adapter to be trimmed" TYPE STRING
# RUNTIME R-4.1.1-asv

# ES 13.10.2022
# PARAMETER paired: "Is the data paired end or single end reads" TYPE [paired, single] DEFAULT single (Are all the reads paired end, so one forward and one reverse FASTQ file for one sample. If single end reads,use only those forward parameters.)
# INPUT OPTIONAL reads.tar: "Tar package containing the FASTQ files" TYPE GENERIC

source(file.path(chipster.common.path,"tool-utils.R"))
source(file.path(chipster.common.path,"zip-utils.R"))

# load 2 libraries
library(Biostrings)
library(ShortRead)

if (adapter3=="" && adapter5==""){
  stop(paste('CHIPSTER-NOTE: ',"You need to give at least one adapter/primer sequence as a parameter"))
}

#check out if the file is compressed and if so unzip it
unzipIfGZipFile("file1.fastq")

#check out if the file is compressed and if so unzip it
#unzipIfGZipFile("reads.tar")

# Make input  folders 
#system("mkdir input_folder")

# untar the tar package to input_folder and list the filenames
#untar("reads.tar", exdir = "input_folder")

#list the full file names with the relative path
#filenames <- sort(list.files("input_folder", full.names=TRUE))

#forward_file=filenames[0]
#if (paired=="paired"){
#    forward_file=filenames[1]
#}

# to ensure we have the right primers, and the correct orientation of the primers 
allOrients <- function(primer) {
    # Create all orientations of the input sequence
    require(Biostrings)
    dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
    orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
        RevComp = reverseComplement(dna))
    return(sapply(orients, toString))  # Convert back to character vector
}
# make a summary txt file for all the different orientations. User can copy these and use in next step
sink("orientations_summary.txt")
    if (!is.na(adapter5)){
        cat("All orientations of the 5' adapter\n")
        orients5 <-allOrients(adapter5)
        print(orients5)
    }
    if (!is.na(adapter5)){
        cat("\nAll orientations of the 3' adapter:\n")
        orients3 <-allOrients(adapter3)
        print(orients3)
    }
sink()

primerHits <- function(primer, fn) {
    print(fn)
    # Counts number of reads in which the primer is found
    nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
    return(sum(nhits > 0))
}


# make a tsv table how many hits were found and in which orientation
if (fileOk("file2.fastq")){ #paired reads, reverse and forward read to check
    if (adapter5==""){
        matrix <-rbind(Forward_reads.3adapter = sapply(orients3, primerHits,fn= "file1.fastq"),
                Reverse_reads.3adapter = sapply(orients3, primerHits, fn = "file2.fastq"))
    }else if (adapter3==""){
        matrix <-rbind(Forward_reads.5adapter = sapply(orients5, primerHits, fn = "file1.fastq"),
            Reverse_reads.5adapter = sapply(orients5, primerHits, fn = "file2.fastq"))
    }else{
        matrix <-rbind(Forward_reads_5adapter = sapply(orients5, primerHits, fn = "file1.fastq"),
                Forward_reads_3adapter = sapply(orients3, primerHits, fn = "file1.fastq"),
                Reverse_reads_5adapter = sapply(orients5, primerHits, fn = "file2.fastq"),
                Reverse_reads_3adapter = sapply(orients3, primerHits, fn = "file2.fastq"))
    }
}else{ #single reads so no reverse read here to specify
    if (adapter5==""){
        matrix <- rbind(Forwardreads_3adapter = sapply(orients3, primerHits, fn = "file1.fastq"))
    }else if (adapter3==""){
        matrix <- rbind(Forward_reads_5adapter = sapply(orients5, primerHits, fn = "file1.fastq"))
    }else{
    
        matrix <- rbind(Forward_reads_5adapter = sapply(orients5, primerHits, fn = "file1.fastq"),
                        Forward_reads_3adapter = sapply(orients3, primerHits, fn = "file1.fastq"))
        print("JEE")
    }
}

write.table(matrix, file ="primer_summary.tsv", sep='\t')

#EOF