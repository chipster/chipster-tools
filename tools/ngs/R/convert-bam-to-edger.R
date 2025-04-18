# TOOL convert-bam-to-edger.R: "Convert genomic BAM file to count table" (This tool takes BAM files as an input, calculates the number of times each sequence tag is identified and removes the ones for which the count is under the user defined threshold.)
# INPUT bam_file.bam: "Alignment against genome in BAM format" TYPE GENERIC
# OUTPUT counts.tsv: "A converted BAM file suitable for DEG analysis"
# PARAMETER count_limit: "Count limit" TYPE INTEGER FROM 0 TO 100000 DEFAULT 0 (The lowest number of times a sequence tag has to appear in the data)
# PARAMETER merge_overlapping: "Merge overlapping reads" TYPE [no, start, end, both] DEFAULT both (If enabled, reads with identical start or end position, but differing in length or nucleotide composition, will be merged into one. The count number will be summed over the merged reads and the sequence and length of the longest among the merged reads will be reported.)

# MG 15.6.2011
# MG, 19.8.2011, added parameter to merge reads with same start or end position
# MK, 26.08.2013, dropping of reads was done at a wrong place
# MK, 10.04.2014, fix bug preveting the use of count_limit filtering for cases where merge_overlapping has been set to no
# MK, 10.04.2014, speeding up the script subtantially by removing for loops

# Extract the BAM file into SAM
samtools.binary <- c(file.path(chipster.tools.path, "samtools-0.1.19", "samtools"))
samtools.command <- paste(samtools.binary, "view bam_file.bam > sam_file")
system(samtools.command)

# Extract the following info from the SAM file:
# nucleotide sequence
# chromosome
# start
# length of sequence
input.file <- "sam_file"
output.file <- "sam_file_extracted"
extract.command <- paste("awk '{print $10\"\t\"$3\"\t\"$4\"\t\"length($10)+$4-1\"\t\"length($10)}'", input.file, ">", output.file)
system(extract.command)

# Sort sequence reads according to chromosome and start position
input.file <- "sam_file_extracted"
output.file <- "sam_file_sorted"
sort.command <- paste("sort -k3n -k4n", input.file, ">", output.file)
system(sort.command)

# Find the unique reads
input.file <- output.file
output.file <- "sam_file_unique"
unique.command <- paste("uniq -c -f 1", input.file, ">", output.file)
system(unique.command)

# Create an output file
input.file <- output.file
output.file <- "sam_file_output"
output.command <- paste("awk '{ print $2\"\t\"$3\"\t\"$4\"\t\"$5\"\t\"$6\"\t\"$1 }' ", input.file, ">", output.file)
system(output.command)

# Creat sequence read ID composed of chromosome name, start position and end position
input.file <- output.file
output.file <- "sam_file_id"
id.command <- paste("awk '{print $2\"_\"$3\"_\"$4\"_\"$1\"\t\"$1\"\t\"$2\"\t\"$3\"\t\"$4\"\t\"$5\"\t\"$6}'", input.file, ">", output.file)
system(id.command)

# Remove sequence reads that were not mapped
input.file <- output.file
output.file <- "sam_file_trimmed"
trim.command <- paste("grep -v \\*", input.file, ">", output.file)
system(trim.command)

# Add column headers
headers <- paste("id\t", "sequence\t", "chr\t", "start\t", "end\t", "length\t", "count", sep = "")
input.file <- output.file
header.file <- "header_file"
output.file <- "counts.tsv"
system(paste("echo \"", headers, "\"", ">", header.file))
merge.command <- paste("cat", header.file, input.file, ">", output.file)
system(merge.command)

# Make final output
file <- c("counts.tsv")
dat <- read.table(file, header = T, sep = "\t", row.names = 1)

if (merge_overlapping == "no") {
    # Don't merge any sequences
    results_table <- dat
}
if (merge_overlapping == "start" || merge_overlapping == "both") {
    # Merge sequences with same starting position into single one and sum counts
    # Go through data chromosome by chromosome
    chr_list <- levels(dat$chr)
    chr_number <- length(chr_list)
    results_table <- NULL

    for (count in 1:chr_number) {
        dat_temp <- dat[dat$chr == chr_list[count], ]
        dat_temp <- droplevels(dat_temp)

        counts_temp <- dat_temp$count
        start_temp <- dat_temp$start
        end_temp <- dat_temp$end
        length_temp <- dat_temp$length
        names(counts_temp) <- rownames(dat_temp)

        counts_temp_merged <- aggregate(counts_temp, by = (list(as.character(dat_temp$start))), FUN = sum)[, 2]
        start_temp_longest <- aggregate(start_temp, by = (list(as.character(dat_temp$start))), FUN = max)[, 2]
        end_temp_longest <- aggregate(end_temp, by = (list(as.character(dat_temp$start))), FUN = max)[, 2]
        length_temp_longest <- aggregate(length_temp, by = (list(as.character(dat_temp$start))), FUN = max)[, 2]

        id_temp <- paste(chr_list[count], "_", start_temp_longest, "_", end_temp_longest, sep = "")
        id_temp_mat <- data.frame(seq = names(counts_temp), id = gsub("[_A-Z]*$", "", names(counts_temp)))
        sequence_temp_matched <- as.character(dat_temp[match(id_temp, id_temp_mat$id), "sequence"])
        sequence_id <- as.character(rownames(dat_temp[match(id_temp, id_temp_mat$id), ]))

        # construct the result table for the chromosome
        results_temp <- data.frame(
            sequence = sequence_temp_matched,
            chr = chr_list[count],
            start = start_temp_longest,
            end = end_temp_longest,
            length = length_temp_longest,
            count = counts_temp_merged,
            row.names = sequence_id
        )
        results_table <- rbind(results_table, results_temp)
    }
    results_table <- results_table[-1, ]
}

if (merge_overlapping == "end" || merge_overlapping == "both") {
    # Merge sequences with same ending position into single one and sum counts
    # Go through data chromosome by chromosome
    if (merge_overlapping == "both") {
        dat <- results_table
    }
    chr_list <- levels(dat$chr)
    chr_number <- length(chr_list)
    results_table <- dat[1, ]
    results_table[1, ] <- 0
    for (count in 1:chr_number) {
        dat_temp <- dat[dat$chr == chr_list[count], ]
        dat_temp <- droplevels(dat_temp)

        counts_temp <- dat_temp$count
        start_temp <- dat_temp$start
        end_temp <- dat_temp$end
        length_temp <- dat_temp$length
        names(counts_temp) <- rownames(dat_temp)

        counts_temp_merged <- aggregate(counts_temp, by = (list(as.character(dat_temp$end))), FUN = sum)[, 2]
        start_temp_longest <- aggregate(start_temp, by = (list(as.character(dat_temp$end))), FUN = max)[, 2]
        end_temp_longest <- aggregate(end_temp, by = (list(as.character(dat_temp$end))), FUN = max)[, 2]
        length_temp_longest <- aggregate(length_temp, by = (list(as.character(dat_temp$end))), FUN = max)[, 2]

        id_temp <- paste(chr_list[count], "_", start_temp_longest, "_", end_temp_longest, sep = "")
        id_temp_mat <- data.frame(seq = names(counts_temp), id = gsub("[_A-Z]*$", "", names(counts_temp)))
        sequence_temp_matched <- as.character(dat_temp[match(id_temp, id_temp_mat$id), "sequence"])
        sequence_id <- as.character(rownames(dat_temp[match(id_temp, id_temp_mat$id), ]))

        # fetch the corresponding seq_id and sequences
        # id_temp <- paste(chr_list[count], "_", start_temp_longest, "_", end_temp_longest, sep="")
        # sequence_id <- character(length(id_temp))
        # for (count_2 in 1:length(id_temp)) {
        # 	sequence_id[count_2] <- names (counts_temp) [grep (id_temp[count_2], names(counts_temp))]
        # }
        # sequence_temp_matched <- dat_temp[sequence_id,]
        # sequence_temp_matched <- as.character(sequence_temp_matched$sequence)

        # construct the result table for the chromosome
        results_temp <- data.frame(
            sequence = sequence_temp_matched,
            chr = chr_list[count],
            start = start_temp_longest,
            end = end_temp_longest,
            length = length_temp_longest,
            count = counts_temp_merged,
            row.names = sequence_id
        )
        results_table <- rbind(results_table, results_temp)
    }
    results_table <- results_table[-1, ]
}

results_table <- results_table[results_table$count >= count_limit, ]
# Writing data to disk
write.table(data.frame(id = rownames(results_table), results_table), file = "counts.tsv", col.names = T, quote = F, sep = "\t", row.names = F)

# EOF
