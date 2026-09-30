# TOOL general_file_processing.R: "Modify text" (This tool can be used to modify txt, tsv, BED and GTF files. For example, you can replace text, or extract rows containing a given text. Note that the search string is interpreted as an extended regular expression, where characters such as . * + ? ^ $ | and brackets have special functions. To match one of them literally, put a backslash in front of it.)
# INPUT input: "Input file" TYPE GENERIC
# OUTPUT OPTIONAL selected.tsv
# OUTPUT OPTIONAL selected.txt
# OUTPUT OPTIONAL selected.bed
# OUTPUT OPTIONAL selected.gtf
# OUTPUT OPTIONAL file_operation.log
# PARAMETER operation: "Operation" TYPE [select: "Select rows with a regular expression", exclude: "Exclude rows with a regular expression", replace: "Replace text", pick_rows: "Select a set of rows from the file" ] DEFAULT replace (Operation to be performed for the selected text or table file)
# PARAMETER OPTIONAL sstring: "Search string" TYPE UNCHECKED_STRING (Search expression)
# PARAMETER OPTIONAL rstring: "Replacement string" TYPE STRING (Replacement string)
# PARAMETER OPTIONAL startrow: "First row to select" TYPE INTEGER DEFAULT 1 (Number of the first row to be selected. Note that in table files, the header row is considered as the first row.)
# PARAMETER OPTIONAL stoprow: "Last row to select" TYPE INTEGER DEFAULT 10000000 (Number of the last row to be selected.)
# PARAMETER OPTIONAL fstyle: "Input file format" TYPE [txt: Text, tsv: Table, bed: BED, gtf: GTF] DEFAULT txt (Is the input file a text file, tab-delimited table, BED file or GTF file)
# PARAMETER OPTIONAL save_log: "Collect a log file" TYPE [yes: Yes, no: No] DEFAULT no (Collect a log file about the analysis run.)
# RUNTIME R-4.5.1
# TOOLS_BIN ""

# KM 10.4.2015
# AMS 9.11.2015 Added support for compressed input files
# EK 30.8.2023 Possibility to force match to the beginning of a string
# TH 26.9.2026 Rewrote in plain R. The search string is an UNCHECKED_STRING and
#              must never be pasted into a system() command.

chunk.size <- 100000

# read the next chunk of lines, stopping the job if the file seems to be binary
readChunk <- function(con) {
    withCallingHandlers(
        readLines(con, n = chunk.size),
        warning = function(w) {
            if (grepl("embedded nul", conditionMessage(w))) {
                # call. = FALSE, so that the message starts with "Error: " like
                # a top-level stop() and the comp finds it
                stop(paste("CHIPSTER-NOTE:", "The input file seems to be binary, not text"), call. = FALSE)
            }
            # e.g. "incomplete final line"
            invokeRestart("muffleWarning")
        }
    )
}

if (nchar(sstring) > 50) {
    stop(paste("CHIPSTER-NOTE:", "Too long search string"))
}

if (nchar(rstring) > 50) {
    stop(paste("CHIPSTER-NOTE:", "Too long replacement string"))
}

is.text.operation <- operation %in% c("select", "exclude", "replace")

if (is.text.operation) {
    if (nchar(sstring) == 0) {
        stop(paste("CHIPSTER-NOTE:", "Please give a search string"))
    }
    # fail early with a readable message, not in the middle of the file
    tryCatch(
        suppressWarnings(grepl(sstring, "")),
        error = function(e) {
            stop(paste("CHIPSTER-NOTE:", "Invalid regular expression:", sstring, "-", conditionMessage(e)), call. = FALSE)
        }
    )
}

# the comp writes the script in UTF-8, but R may read it in the C locale
Encoding(rstring) <- "UTF-8"

if (is.na(startrow)) {
    startrow <- 1
}
if (is.na(stoprow)) {
    stoprow <- Inf
}

# file() decompresses gzip input on its own when opened in text mode

# Match characters in UTF-8 text, but fall back to bytes if the input isn't
# valid UTF-8, because R would fail on it otherwise. Decide this once for the
# whole file, so that the same expression means the same thing on every row.
use.bytes <- FALSE
if (is.text.operation) {
    input.con <- file("input", open = "r")
    repeat {
        lines <- readChunk(input.con)
        if (length(lines) == 0) {
            break
        }
        if (!all(validUTF8(lines))) {
            use.bytes <- TRUE
            break
        }
    }
    close(input.con)
}

# process the file in chunks to keep the memory usage low also for large GTF files
input.con <- file("input", open = "r")
output.con <- file(paste("selected.", fstyle, sep = ""), open = "w")
rows.in <- 0
rows.out <- 0

repeat {
    lines <- readChunk(input.con)
    if (length(lines) == 0) {
        break
    }
    row.numbers <- rows.in + seq_along(lines)
    rows.in <- rows.in + length(lines)

    if (is.text.operation && !use.bytes) {
        Encoding(lines) <- "UTF-8"
    }

    if (operation == "select") {
        lines <- lines[grepl(sstring, lines, useBytes = use.bytes)]
    } else if (operation == "exclude") {
        lines <- lines[!grepl(sstring, lines, useBytes = use.bytes)]
    } else if (operation == "replace") {
        lines <- gsub(sstring, rstring, lines, useBytes = use.bytes)
    } else if (operation == "pick_rows") {
        lines <- lines[row.numbers >= startrow & row.numbers <= stoprow]
    }

    writeLines(lines, output.con, useBytes = TRUE)
    rows.out <- rows.out + length(lines)

    # the log reports the number of rows in the whole file, so read it all in that case
    if (operation == "pick_rows" && rows.in >= stoprow && save_log == "no") {
        break
    }
}

close(input.con)
close(output.con)

if (save_log == "yes") {
    log.lines <- paste("Operation:", operation)
    if (is.text.operation) {
        log.lines <- c(log.lines, paste("Search string:", sstring))
    }
    if (operation == "replace") {
        log.lines <- c(log.lines, paste("Replacement string:", rstring))
    }
    if (operation == "pick_rows") {
        last.row <- if (is.finite(stoprow)) format(stoprow, scientific = FALSE) else "end of file"
        log.lines <- c(log.lines, paste("Rows:", format(startrow, scientific = FALSE), "-", last.row))
    }
    log.lines <- c(
        log.lines,
        paste("Number of rows in the input file:", format(rows.in, scientific = FALSE)),
        paste("Number of rows in the output file:", format(rows.out, scientific = FALSE))
    )
    writeLines(log.lines, "file_operation.log", useBytes = TRUE)
}
