# TOOL ncbi_blastn.R: "Nucleotide BLAST" (NCBI BLAST for searching a nucleotide database using a nucleotide query. Your query sequence set can contain 10 sequeces in maximum.)
# INPUT query.fa: "Query sequences" TYPE GENERIC
# OUTPUT OPTIONAL blast_results.txt
# OUTPUT OPTIONAL blast_results.xml
# OUTPUT OPTIONAL blast_results.tsv
# OUTPUT OPTIONAL blast_results.fasta
# OUTPUT OPTIONAL blast_results.csv
# OUTPUT OPTIONAL blast_results.asn1
# OUTPUT OPTIONAL blast_results.html
# OUTPUT OPTIONAL blast.log
# PARAMETER db: "Database" TYPE [nt: "NCBI nucleotide collection (nr/nt)", refseq_rna: "Reference RNA sequences (refseq_rna)",  refseq_genomic: "Reference genomic sequences (refseq_genomic)", pdb: "PDB sequences"] DEFAULT nt (Database to search.)
# PARAMETER task: "BLAST program to use" TYPE [blastn: "blastn", blastn-short: "blastn-short",  dc-megablast: "discontiguous megablast",  megablast: "megablast", rmblastn: "rmblastn"] DEFAULT megablast (BLAST algorithm to use. Megablast is intended for comparing a query to closely related sequences and works best if the target percent identity is 95% or more but is very fast. Discontiguous megablast uses an initial seed that ignores some bases (allowing mismatches\) and is intended for cross-species comparisons. BlastN is slow, but allows a word-size down to seven bases.)
# PARAMETER OPTIONAL evalue: "Expectation threshold for saving hits" TYPE DECIMAL DEFAULT 10 (E-value specifies the statistical significance threshold for reporting matches against database sequences. The default value 10 means that 10 such matches are expected to be found merely by chance. Lower thresholds are more stringent, leading to fewer chance matches being reported.)
# PARAMETER OPTIONAL word_size: "Word size" TYPE INTEGER DEFAULT 28 (The length of the seed that initiates an alignment. BLAST works by finding word-matches between the query and database sequences. One may think of this process as finding hot-spots that BLAST can then use to initiate extensions that might eventually lead to full-blown alignments. For nucleotide-nucleotide searches an exact match of the entire word is required before an extension is initiated, so that one normally regulates the sensitivity and speed of the search by increasing or decreasing the word-size.)
# PARAMETER OPTIONAL num_hits: "Maximun number of hits to collect per sequence" TYPE INTEGER DEFAULT 100 (Number of database sequences to show one-line descriptions for.)
# PARAMETER OPTIONAL outfmt: "Output format type" TYPE [0: "normal BLAST report with pairwise alignments", 1: "query-anchored alignments showing identities", 2: "query-anchored alignments with no identities", 3: "flat query-anchored, show identities", 4: "flat query-anchored, no identities", 5: "XML Blast output", 6: "tabular", 10: "comma-separated values", 11: "BLAST archive format", 12: "list of unique hit sequence identifiers", 14: "hit regions in fasta format"] DEFAULT 0 (Output format type)
# PARAMETER OPTIONAL dust: "Filter low complexity regions" TYPE [yes: yes, no: no] DEFAULT yes (Use the DUST program for filtering low complexity regions in the query sequence.)
# PARAMETER OPTIONAL entrez_query: "Entrez query to limit search" TYPE STRING DEFAULT "none" (You can use Entrez query syntax to search a subset of the selected BLAST database. This can be helpful to limit searches to molecule types, sequence lengths or to exclude organisms.)
# PARAMETER OPTIONAL query_loc: "Location on the query sequence" TYPE STRING DEFAULT "full length" (Location of the search region in the query sequence, for example: 23-66.)
# PARAMETER OPTIONAL reward: "Reward for a match" TYPE STRING DEFAULT "default" (Reward for a nucleotide match. The scoring system consists of a reward for a match and a penalty for a mismatch. The absolute reward/penalty ratio should be increased as one looks at more divergent sequences. A ratio of 0.33 (1/-3\) is appropriate for sequences that are about 99% conserved; a ratio of 0.5 (1/-2\) is best for sequences that are 95% conserved; a ratio of about one (1/-1\) is best for sequences that are 75% conserved.)
# PARAMETER OPTIONAL penalty: "Penalty for a mismatch" TYPE STRING DEFAULT "default" (Penalty for a nucleotide mismatch. The scoring system consists of a reward for a match and a penalty for a mismatch. The absolute reward/penalty ratio should be increased as one looks at more divergent sequences. A ratio of 0.33 (1/-3\) is appropriate for sequences that are about 99% conserved; a ratio of 0.5 (1/-2\) is best for sequences that are 95% conserved; a ratio of about one (1/-1\) is best for sequences that are 75% conserved.)
# PARAMETER OPTIONAL gapopen: "Gap opening penalty" TYPE STRING DEFAULT "default" (Gap opening penalty ranging from 6 to 25. Note that if you assign this value, you must define also the gap extension penalty.)
# PARAMETER OPTIONAL gapextend: "Gap extension penalty" TYPE STRING DEFAULT "default" (Gap extension penalty ranging from 1 to 3. Note that if you assign this value, you must define also the gap opening penalty.)
# PARAMETER OPTIONAL save_log: "Output a log file" TYPE [yes: yes, no: no] DEFAULT no (Collect a log file for the BLAST run.)

# KM 31.10.2013
# EK 18.11.2013 changes to parameter names and descriptions

# check out if the file is compressed and if so unzip it
source(file.path(chipster.common.lib.path, "zip-utils.R"))
unzipIfGZipFile("query.fa")
# unzipIfGZipFile("dbprot.fa")

# pb settings
pb.binary <- file.path(chipster.module.path, "/shell/pb_for_chipster.sh")
# pb.binary <- file.path(chipster.tools.path, "blast", "/ncbi-blast-2.2.29+", "bin", "pb_for_chipster")
command.start <- paste(pb.binary, "blastn")

emboss.path <- file.path(chipster.tools.path, "emboss", "bin")



# check sequece file type
inputfile.to.check <- ("query.fa")
sfcheck.binary <- file.path(chipster.module.path, "/shell/sfcheck.sh")
sfcheck.command <- paste(sfcheck.binary, emboss.path, inputfile.to.check)
str.filetype <- system(sfcheck.command, intern = TRUE)

if (str.filetype == "Not an EMBOSS compatible sequence file") {
    stop("CHIPSTER-NOTE: Your input file is not a sequence file that is compatible with the tool you try to use")
}


# count the query sequeces
seqcount.exe <- file.path(emboss.path, "seqcount -filter query.fa")
str.queryseq <- system(seqcount.exe, intern = TRUE)
num.queryseq <- as.integer(str.queryseq)


# In blast command outputformat is used as a string
# in job control as a number
outfmt.string <- outfmt
outfmt <- as.integer(outfmt)
# round(num.queryseq)

if (num.queryseq > 10) {
    stop(paste("CHIPSTER-NOTE: Too many query sequences. Maximun is 10 but your file contains ", num.queryseq))
}

# Modify table format
outfmt.is.table <- paste("no")

if (outfmt == 6) {
    outfmt.string <- paste('"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle"')
    outfmt.is.table <- paste("yes")
}
# These parameters have allways some value
general.parameters <- paste("-chipster_path /opt/chipster -remote -no_slurm -query query.fa -out blast_results -db", db)
general.parameters <- paste(general.parameters, "-task", task)
general.parameters <- paste(general.parameters, "-evalue ", evalue)
general.parameters <- paste(general.parameters, "-dust", dust)
general.parameters <- paste(general.parameters, "-outfmt", outfmt.string)
general.parameters <- paste(general.parameters, "-num_threads", chipster.threads.max)

optional.parameters <- paste(" ")

if (nchar(entrez_query) > 4) {
    optional.parameters <- paste(optional.parameters, '-entrez_query "', entrez_query, '"')
}


if (outfmt < 5) {
    optional.parameters <- paste(optional.parameters, " -html -num_descriptions ", num_hits, "-num_alignments", num_hits)
}

if (outfmt > 5) {
    optional.parameters <- paste(optional.parameters, " -max_target_seqs ", num_hits)
}

# Check text formatted parameters
if (nchar(gapopen) > 2) {
    gapopen <- paste("default")
}
if (nchar(gapopen) < 1) {
    gapopen <- paste("default")
}

if (nchar(gapextend) > 2) {
    gapextend <- paste("default")
}
if (nchar(gapextend) < 1) {
    gapextend <- paste("default")
}

if (nchar(reward) < 1) {
    reward <- paste("default")
}

if (nchar(reward) > 2) {
    reward <- paste("default")
}

if (nchar(penalty) < 1) {
    penalty <- paste("default")
}

if (nchar(penalty) > 2) {
    penalty <- paste("default")
}

if (nchar(query_loc) < 3) {
    penalty <- paste("full length")
}

# set up more parameters if applied

if (gapopen != "default") {
    if (gapextend == "default") {
        stop(paste("If you define gap opening penalty you must also define gap extension pelalty"))
    }
}

if (gapextend != "default") {
    if (gapopen == "default") {
        stop(paste("If you define gap opening penalty you must also define gap extension pelalty"))
    }
}

if (gapopen != "default") {
    if (gapextend != "default") {
        optional.parameters <- paste(optional.parameters, " -gapopen ", gapopen, "-gapextend", gapextend)
    }
}


if (query_loc != "full length") {
    optional.parameters <- paste(optional.parameters, " -query_loc ", query_loc)
}

if (reward != "default") {
    optional.parameters <- paste(optional.parameters, " -reward ", reward)
}

if (penalty != "default") {
    optional.parameters <- paste(optional.parameters, " -penalty ", penalty)
}
command.end <- (paste(" >>", "blast.log 2>&1"))

# run the task
system("ls -l > blast.log")
blast.command <- paste(command.start, general.parameters, optional.parameters, command.end)
echo.command <- paste("echo '", blast.command, "' >> blast.log")
system(echo.command)
system(blast.command)
system("echo xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx")
system("ls -l >> blast.log")
# rename the results
if (outfmt < 5) {
    test.command <- paste("echo renaming as html. outformat is ", outfmt, ">> blast.log")
    system(test.command)
    system("mv blast_results blast_results.html")
}

if (outfmt == 5) {
    system("mv blast_results blast_results.xml")
}

if (outfmt.is.table == "yes") {
    system('printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "query id" "subject id" "identity percent" "alignment length" mismatches "gap opens" "query start" "query end" "subject start" "subject end" evalue "bit score" "subject description"> blast_results.tsv ')
    system("cat blast_results >> blast_results.tsv")
}

if (outfmt == 7) {
    system("mv blast_results blast_results.tsv")
}

if (outfmt == 8) {
    system("mv blast_results blast_results.txt")
}

if (outfmt == 10) {
    system("mv blast_results blast_results.csv")
}

if (outfmt == 11) {
    system("mv blast_results blast_results.asn1")
}

if (outfmt == 12) {
    system("mv blast_results blast_results.txt")
}

if (outfmt == 13) {
    system("mv blast_results blast_results.fasta")
}

if (outfmt == 14) {
    system("mv blast_results blast_results.fasta")
}

if (save_log == "no") {
    system("rm -f blast.log")
}
