# TOOL dimont-extract-data.R: "Dimont sequence extractor" (Extracts genomic regions specified in a BED-like file format in the annotated FastA format as required by Dimont.)
# INPUT regions.bed: "Genomic regions" TYPE GENERIC (The genomic regions to be extracted in a BED-like file format, e.g., BED, GTF, narrowPeak.)
# OUTPUT extracted.fasta: "Extracted sequences" (The sequences extracted from the given genome using the supplied region specifications.)
# PARAMETER organism: "Genome" TYPE ["FILES genomes/fasta .fa"] DEFAULT "SYMLINK_TARGET genomes/fasta/default .fa" (Genome that the regions refer to.)
# PARAMETER chromcol: "Chromosome column" TYPE INTEGER FROM 1 DEFAULT 1 (The column of the Regions file, which contains the chromosome information.)
# PARAMETER startcol: "Start column" TYPE INTEGER FROM 1 DEFAULT 2 (The column of the Regions file containing the start position of the genomic region.)
# PARAMETER seccoord: "Second coordinate" TYPE INTEGER FROM 1 DEFAULT 3 (The second genomic coordinate with meaning specified by parameter \"Meaning of second coordinate\".)
# PARAMETER seccol: "Meaning of second coordinate" TYPE [center: "Center of peak (relative to start\)", end: "End of peak (global coordinates\)"] DEFAULT end (The meaning of the second genomic coordinate. This may either be the position of the peak summit relative to the position in Start, or the end position of the peak.)
# PARAMETER statcol: "Statistics column" TYPE INTEGER FROM 1 DEFAULT 5 (The column containing the peak statistics information (or another measure of peak confidence\).)
# PARAMETER OPTIONAL width: Width TYPE INTEGER FROM 1 DEFAULT 1000 (The width of the genomic region to be extracted. Recommended values: 1000 for ChIP-seq and 100 for ChIP-exo.)
# PARAMETER OPTIONAL has.row.names: "Regions file has row names" TYPE [yes, no] DEFAULT no (BED file has row names in its first column.)


# MK: 11.12.2013. Since also row-names are printed, column indexes had to be shifted by one
# JG: 30.01.2014: Code for chr=="no" does not work, data extraction script works in either case. Hence, parameter and code removed
# EK: 3.6.2014: removed chr parameter that was commented out already
# AMS 04.07.2014 New genome/gtf/index locations & names

genome.fa <- file.path(chipster.tools.path, "genomes", "fasta", paste(organism, ".fa", sep = "", collapse = ""))

tool <- file.path(chipster.tools.path, "dimont", "extract_data_single_chipster.pl")
if (has.row.names == "yes") {
    command <- paste("perl", tool, genome.fa, "regions.bed", (chromcol + 1), (startcol + 1), seccol, (seccoord + 1), width, (statcol + 1), "extracted.fasta")
} else {
    command <- paste("perl", tool, genome.fa, "regions.bed", (chromcol + 0), (startcol + 0), seccol, (seccoord + 0), width, (statcol + 0), "extracted.fasta")
}

system(command)
