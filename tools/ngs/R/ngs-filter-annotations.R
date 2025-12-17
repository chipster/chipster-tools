# TOOL ngs-filter-annotations.R: "Filter table by column term" (Allows the user to filter the table rows on the basis of terms in any text column.)
# INPUT annotations.tsv: "Table to be filtered"  TYPE GENERIC
# INPUT OPTIONAL genelist.tsv: "List of terms" TYPE GENERIC (List of terms \(like gene names\) to look for as a simple tsv file \(no header, just list of terms, one term in each row\). Note that giving this file as an input overwrites the \"term to match\" parameter!)
# OUTPUT filtered-NGS-results.tsv: filtered-NGS-results.tsv
# PARAMETER column: "Column to filter by" TYPE COLUMN_SEL (Data column to filter by)
# PARAMETER match.term: "Term to match" TYPE STRING DEFAULT empty (Textual term to search for. If you list multiple terms, separate them with comma \(,\). Note that this list is NOT used if you give a list as a tsv file as input!)
# PARAMETER has.rownames: "Does the first column lack a title" TYPE [no: no, yes: yes] DEFAULT no (Specifies whether the first column has a title or not.)
# PARAMETER mode: Mode TYPE [include: include, exclude: exclude] DEFAULT include (Defines whether the found terms\/genes should be included or excluded from the resulting data table.)
# RUNTIME R-4.5.1
# TOOLS_BIN ""

# MG 29.5.2010
# MK, EK 21.08.2013 added support for rownames
# ML 20.4.2021 Add option to input multiple gene names with comma and use an input file, and option to include or exclude the wanted terms.

# Loads the normalized data
file <- c("annotations.tsv")

if (has.rownames == "yes") {
    # dat <- read.table(file, header=TRUE, sep="\t", row.names=NULL, quote = "", comment.char="")
    dat <- read.table(file, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE, quote = "", comment.char = "")
}
if (has.rownames == "no") {
    dat <- read.table(file, header = TRUE, sep = "\t", row.names = NULL, quote = "", comment.char = "")
    # dat <- read.table(file, header=TRUE, sep="\t", row.names=1, check.names=FALSE, quote = "", comment.char="")
}

dat <- read.table(file, header = T, sep = "\t", check.names = FALSE, quote = "", comment.char = "")




# Input query terms:
# If a tsv file for terms (for example gene names) is given, use that.
# A simple .tsv file with just the query terms one per row, no headers, is expected.
if (file.exists("genelist.tsv")) {
    genelist <- read.table("genelist.tsv", header = F, sep = "\t", row.names = 1)
    gene.names.clean <- rownames(genelist)
} else {
    # Otherwise, use the list of genes given as a parameter.
    # Handle first:
    # Split from comma:
    gene.names.split <- strsplit(match.term, ",")[[1]]
    # remove whitespace:
    gene.names.clean <- gsub("[[:blank:]]", "", gene.names.split)
}


# Looking at the column name:
# When looking at the untitled first column, column = " " or "" or "untitled".
# Usually. It's a wild west.
# In example cases we had, it was "", even though print(paste("", column, "")) gives "_ _" so it LOOKS to be " ".
if ((column == " ") || (column == "") || (column == "untitled") || is.null(column)) {
    column <- "first_column_withouth_title"
} else {
    # Otherwise, extract the column data to a vector:
    f <- dat[, colnames(dat) %in% column]
    if (length(which(colnames(dat) %in% column == TRUE)) > 1) {
        stop("CHIPSTER-NOTE: Please choose a column that occurs only once in the table")
    }
}

# If user wants to INCLUDE the rows with the query terms from data:
if (mode == "include") {
    # Loop trough the list of query terms (gene names) if multiple were given:
    for (i in 1:length(gene.names.clean)) {
        match.term <- gene.names.clean[i]

        # When looking at the untitled first column:
        if (column == "first_column_withouth_title") {
            dat2 <- dat[grep(match.term, rownames(dat)), ]
        } else {
            # Extract the data from the column in question
            # dat2 <- dat[grep(match.term, dat[,(as.vector(grep(column,names(dat))))]),]
            dat2 <- dat[grep(match.term, f), ]
        }
        # Merge the tables of different query terms:
        if (i == 1) {
            dat3 <- dat2
        } else {
            dat3 <- rbind(dat3, dat2)
        }
    }
    # Remove duplicate rows (in case some query words find the same rows)
    dat3 <- dat3[!duplicated(dat3), ]
}


# If user wants to EXCLUDE the rows with the query terms from data:
if (mode == "exclude") {
    # Loop trough the list of query terms (gene names) if multiple were given:
    for (i in 1:length(gene.names.clean)) {
        match.term <- gene.names.clean[i]

        # When looking at the untitled first column:
        if (column == "first_column_withouth_title") {
            dat2 <- dat[grep(match.term, rownames(dat), invert = TRUE), ]
        } else {
            # Extract the data from the column in question
            # dat2 <- dat[grep(match.term, dat[,(as.vector(grep(column,names(dat))))], invert=TRUE),]
            dat2 <- dat[grep(match.term, f, invert = TRUE), ]
        }
        # After first round, start from the previously extracted table:
        dat <- dat2
        f <- dat[, colnames(dat) %in% column]
    }
    dat3 <- dat2
}


# Write the data to disk
if (has.rownames == "yes") {
    write.table(dat3, "filtered-NGS-results.tsv", sep = "\t", row.names = T, col.names = T, quote = F)
}
if (has.rownames == "no") {
    write.table(dat3, "filtered-NGS-results.tsv", sep = "\t", row.names = F, col.names = T, quote = F)
}
