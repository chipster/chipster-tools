read_input_definitions <- function() {
    # Read in the data
    list <- scan("chipster-inputs.tsv", what = "", sep = "\n", comment.char = "#")
    # Separate elements by one or more whitepace
    inputdef <- strsplit(list, "[\t]+")
    # Extract the first vector element and set it as the list element name
    names(inputdef) <- sapply(inputdef, function(list) list[[1]])
    # Remove the first vector element from each list element
    inputdef <- lapply(inputdef, function(list) list[-1])

    return(inputdef)
}

# Reads chipster-inputs.tsv as a data frame with the input name in column 1 and the dataset name in
# column 2, both as character: dataset names can look like numbers, or be "NA".
#
read_input_names <- function() {
    # comp writes comment lines before the names (ToolUtils.writeInputDescription). They are skipped by
    # count: a dataset name can contain #, so comment.char can't be used. read.table(text = ...) would
    # change non-ASCII bytes to <c3><a4> in the C locale, so the file is read again.
    is.comment <- grepl("^#", readLines("chipster-inputs.tsv", warn = FALSE), useBytes = TRUE)
    skip <- sum(cumprod(is.comment))
    read.table(
        "chipster-inputs.tsv",
        header = FALSE,
        sep = "\t",
        skip = skip,
        colClasses = "character",
        na.strings = character(0),
        quote = "",
        comment.char = ""
    )
}

write_output_definitions <- function(output_names) {
    write.table(
        output_names,
        file = "chipster-outputs.tsv",
        row.names = FALSE,
        col.names = FALSE,
        quote = FALSE,
        sep = "\t"
    )
}

# Removes user-specified postfix if present
#
remove_postfix <- function(name, postfix) {
    pf <- paste(postfix, "$", sep = "")
    if (grepl(pf, name)) {
        basename <- substr(name, 1, (nchar(name) - nchar(postfix)))
        return(basename)
    } else {
        return(name)
    }
}

# Removes extension, i.e. everything after the last dot (including the dot)
#
remove_extension <- function(name) {
    return(sub("\\.[^.]*$", "", name))
}

# Strips common file extensions from a file name
#
strip_name <- function(name) {
    known_postfixes <- c(
        ".gz",
        ".bam",
        ".sam",
        ".fa",
        ".fasta",
        ".fq",
        ".fastq",
        ".gtf",
        ".tsv",
        ".txt",
        "vcf",
        "_trimmed",
        "_filtered"
    )
    newname <- name
    while (TRUE) {
        for (i in known_postfixes) {
            newname <- remove_postfix(newname, i)
        }
        if (nchar(newname) == nchar(name)) {
            break
        }
        name <- newname
    }
    return(name)
}

# Tells which strings are valid UTF-8, like validUTF8(). validUTF8() came in R 3.3.0, but the default
# runtime R is R 3.2.3. There the strings are checked with utf8ToInt(), which gives NA or an error for
# invalid UTF-8. It accepts code points above U+10FFFF, so they are rejected separately. Unlike
# validUTF8(), it rejects the noncharacters U+FFFE and U+FFFF, which are valid UTF-8 but never needed
# in names. NA is valid, as in validUTF8().
#
valid_utf8 <- if (exists("validUTF8", mode = "function", envir = baseenv())) {
    validUTF8
} else {
    function(x) {
        vapply(
            x,
            function(s) {
                if (is.na(s)) {
                    return(TRUE)
                }
                code.points <- tryCatch(utf8ToInt(s), error = function(e) NA)
                !anyNA(code.points) && all(code.points <= 0x10FFFF)
            },
            logical(1),
            USE.NAMES = FALSE
        )
    }
}

# Checks names that come from the user's data, e.g. the file names in a tar package or the
# lines of a list file, before they are used in commands. Allows only letters (also ä, ö, å etc.),
# numbers and . _ + -, so that the shell sees each name as one plain word. Stops the job otherwise.
# Returns the names, so that the check can wrap the listing: names <- safe_name(list.files(...))
#
safe_name <- function(names) {
    # Encoding<- fails on an empty vector
    if (length(names) == 0) {
        return(names)
    }
    # The job may run in the C locale, where R would see the bytes of ä instead of a letter
    utf8 <- names
    Encoding(utf8)[Encoding(utf8) == "unknown"] <- "UTF-8"
    valid <- !is.na(utf8) & valid_utf8(utf8)
    # \p{M} for the accents of decomposed letters, which e.g. macOS uses in file names
    valid[valid] <- grepl("^[\\p{L}\\p{M}\\p{N}._+][\\p{L}\\p{M}\\p{N}._+-]*\\z", utf8[valid], perl = TRUE)
    unsafe <- names[!valid | names %in% c(".", "..")]
    if (length(unsafe) > 0) {
        stop(paste(
            "CHIPSTER-NOTE: File and sample names can contain only letters, numbers and the characters . _ + -",
            "and they can't start with -. Please rename these:",
            paste(unsafe, collapse = ", ")
        ))
    }
    return(names)
}

# Checks names from the user's data that are used only as file names in R or in quoted commands, where
# safe_name() would reject too much. Allows everything except names that would point outside the folder.
# Stops the job otherwise. Returns the names invisibly, so that a plain check call doesn't print them.
#
safe_file_name <- function(names) {
    # R expands a leading ~ to the home folder
    unsafe <- names[
        is.na(names) |
            grepl("/", names, fixed = TRUE, useBytes = TRUE) |
            grepl("^~", names, useBytes = TRUE) |
            names %in% c("", ".", "..")
    ]
    if (length(unsafe) > 0) {
        stop(paste(
            "CHIPSTER-NOTE: Names can't be empty, . or .., start with ~ or contain /.",
            "Please rename these:",
            paste0("\"", unsafe, "\"", collapse = ", ")
        ))
    }
    invisible(names)
}

# Reads the lines of a list file from the user. Removes the byte order mark that e.g. Windows Notepad
# may add to the beginning of the file, so that it doesn't end up in the first name. Checks the bytes,
# because the job may run in the C locale. readLines(filename) reads also gzipped files. Don't use
# gzfile(): in R 3.2.3 it drops the last line if the file doesn't end with a newline. Skips empty
# lines and lines with only whitespace, e.g. a blank line at the end of the file. Stops with a note
# if the file is in UTF-16, which e.g. Notepad calls "Unicode". Tests are in tests/test-tool-utils.R.
#
read_list_file <- function(filename) {
    # warn = FALSE: a list file without a final newline is fine
    lines <- readLines(filename, warn = FALSE)
    if (length(lines) > 0) {
        first <- charToRaw(lines[1])
        # readLines() cuts the lines of a UTF-16 file at its nul bytes, so only the byte order mark is
        # left to recognize it from
        if (
            length(first) >= 2 && (all(first[1:2] == as.raw(c(0xff, 0xfe))) || all(first[1:2] == as.raw(c(0xfe, 0xff))))
        ) {
            stop(paste(
                "CHIPSTER-NOTE: The list file",
                filename,
                "is in UTF-16 encoding. Please save it as UTF-8,",
                "e.g. in Notepad with Save as, Encoding: UTF-8."
            ))
        }
        if (length(first) >= 3 && all(first[1:3] == as.raw(c(0xef, 0xbb, 0xbf)))) {
            lines[1] <- rawToChar(first[-(1:3)])
        }
    }
    # useBytes, because a line may not be valid in the job's locale. \xc2\xa0 is a non-breaking space in
    # UTF-8, which e.g. Word may leave on an empty line.
    lines <- lines[!grepl("^([[:space:]]|\xc2\xa0)*$", lines, useBytes = TRUE)]
    return(lines)
}

# If the names look like typical paired-end names: *_1, *_2, remove the ending and return the name.
# If not, return first name as-is.
#
paired_name <- function(name1, name2) {
    if (grepl("_1$", name1) && grepl("_2$", name2)) {
        return(remove_postfix(name1, "_1"))
    }
    if (grepl("_2$", name1) && grepl("_1$", name2)) {
        return(remove_postfix(name1, "_2"))
    }
    return(name1)
}

# Takes as an argument a file name of a list of input names. It Goes through the list and compares
# it to chipster-inputs.tsv. Returns a vector of internal input names or gives an error message,
# if the list is empty or a listed input file is missing.
#
make_input_list <- function(listfile) {
    # read input names
    input.names <- read_input_names()

    # read list file, without a byte order mark and blank lines
    name.list <- read_list_file(listfile)
    if (length(name.list) == 0) {
        # The user knows the list by its dataset name
        list.name <- input.names[match(listfile, input.names[, 1]), 2]
        if (is.na(list.name)) {
            list.name <- listfile
        }
        stop(paste("CHIPSTER-NOTE: The input file list", list.name, "is empty. List one file name per line."))
    }

    # Check for duplicated etries
    if (anyDuplicated(name.list)) {
        message <- paste(
            "Input file list contains one or more duplicated entries:",
            paste(unique(name.list[duplicated(name.list)]), collapse = ", "),
            "\nFilenames must be unique so they can be correctly assigned."
        )
        stop(paste("CHIPSTER-NOTE: ", message))
    }

    # Check that inputs exist
    sdf <- setdiff(name.list, input.names[, 2])
    if (identical(sdf, character(0))) {
        # A listed name must identify one input
        dataset.names <- input.names[, 2]
        ambiguous <- intersect(name.list, dataset.names[duplicated(dataset.names)])
        if (length(ambiguous) > 0) {
            message <- paste(
                "Several selected files have the same name:",
                paste(ambiguous, collapse = ", "),
                "\nRename them so that they can be told apart in the list file."
            )
            stop(paste("CHIPSTER-NOTE: ", message))
        }
        # Exact match, dataset names can contain regex characters like ( ) + .
        input.list <- input.names[match(name.list, dataset.names), 1]
    } else {
        message <- paste(
            "Input file list includes one or more files that has not been selected:",
            paste(sdf, collapse = ", ")
        )
        stop(paste("CHIPSTER-NOTE: ", message))
    }
    return(input.list)
}

# fileOk and fileNotOk ar functions to check that file exists and optionally that it is larger than the provided threshold.
# fileOk("myfile.txt") : returns TRUE if file myfile.txt exists
# fileOk("myfile.txt", 100) : returns TRUE if file myfile.txt exists and is min. 100 bytes
# fileOk("myfile.txt",,2) : returnd TRUE if file myfile.txt exists and is min. 2 lines
# The two functions work the same exept fileOk returns TRUE if the file meets the criteria, and fileNotOk returns TRUE
# if it does not. The idea is to make if statements more readable.

fileOk <- function(filename, minsize, minlines) {
    # Check that file exists
    if (!(file.exists(filename))) {
        return(FALSE)
    }
    # Check minimum size if provided
    if (!(missing(minsize))) {
        if (file.info(filename)$size < minsize) {
            return(FALSE)
        }
    }
    # Check line number if provided
    if (!(missing(minlines))) {
        linenumber <- as.integer(system(paste("wc -l <", filename), intern = TRUE))
        if (linenumber < minlines) {
            return(FALSE)
        }
    }
    # If all checks pass, return TRUE
    return(TRUE)
}

fileNotOk <- function(filename, minsize, minlines) {
    # Check if file exists
    if (!(file.exists(filename))) {
        return(TRUE)
    }
    # Check minimum size if provided
    if (!(missing(minsize))) {
        if (file.info(filename)$size < minsize) {
            return(TRUE)
        }
    }
    # Checkline number if provided
    if (!(missing(minlines))) {
        linenumber <- as.integer(system(paste("wc -l <", filename), intern = TRUE))
        if (linenumber < minlines) {
            return(TRUE)
        }
    }
    # If all checks pass, return FALSE
    return(FALSE)
}

fileCheck <- function(filename, minsize, minlines) {
    # Check if file exists
    if (!(file.exists(filename))) {
        stop(paste("CHIPSTER-NOTE: ", "Required file", filename, "does not exist."))
    }
    # Check minimum size if provided
    if (!(missing(minsize))) {
        filesize <- file.info(filename)$size
        if (filesize < minsize) {
            stop(paste(
                "CHIPSTER-NOTE: ",
                "Required file",
                filename,
                "exists, but is smaller (",
                filesize,
                ") than required size (",
                minsize,
                ")."
            ))
        }
    }
    # Checkline number if provided
    if (!(missing(minlines))) {
        linenumber <- as.integer(system(paste("wc -l <", filename), intern = TRUE))
        if (linenumber < minlines) {
            stop(paste(
                "CHIPSTER-NOTE: ",
                "Required file",
                filename,
                "exists, but has fewer lines (",
                linenumber,
                ") than required (",
                minlines,
                ")."
            ))
        }
    }
}

# Returns TRUE if file is in FASTA format. File needs to be uncompressed.
isFasta <- function(filename) {
    emboss.path <- file.path(chipster.tools.path, "emboss", "bin")
    sfcheck.binary <- file.path(chipster.module.path, "../misc/shell/sfcheck.sh")
    sfcheck.command <- paste(sfcheck.binary, emboss.path, shQuote(filename))
    str.filetype <- system(sfcheck.command, intern = TRUE)
    if (str.filetype == "fasta") {
        return(TRUE)
    } else {
        return(FALSE)
    }
}

# Returns TRUE if file is in FASTQ format. File needs to be uncompressed.
isFastq <- function(filename, version = 1) {
    # Choose Emboss version according to runtime
    if (version == 1) {
        emboss.path <- file.path(chipster.tools.path, "emboss", "bin")
    } else {
        emboss.path <- file.path(chipster.tools.path, "emboss-24.04", "bin")
    }

    sfcheck.binary <- file.path(chipster.module.path, "../misc/shell/sfcheck.sh")
    sfcheck.command <- paste(sfcheck.binary, emboss.path, shQuote(filename))
    str.filetype <- system(sfcheck.command, intern = TRUE)
    if (grepl("fastq", str.filetype)) {
        return(TRUE)
    } else {
        return(FALSE)
    }
}

# Wrapper for system2() command.
# Optionally captures stderr.
# Checks for exit status and gives an error message if status != 0
# Allows setting environment variables give as character vector in style "VARIABLE=value"
#
runExternal <- function(command, env = NULL, capture = TRUE, checkexit = TRUE) {
    # Split command to words
    wcom <- strsplit(command, " ")[[1]]

    # if environment varaiables are set, the command has to be started with "bash -c" and encased in single quotes.
    if (!is.null(env)) {
        wcom <- c("bash", "-c", "\'", wcom, "\'")
    }

    # Add timestamp and command line to stderr.log
    system("echo \"# Chipster log\" >> stderr.log")
    system("echo \"# \"\`date\` >> stderr.log")
    write(paste(wcom, collapse = " "), file = "stderr.log", append = TRUE)

    # Capture of stderr is optional
    if (capture) {
        # Run command, capture stderr
        exitcode <- system2(wcom[1], wcom[2:length(wcom)], stderr = "stderr.tmp", env = env)

        # Append error messages to the log file. stderr.tmp is overwritten each time runExternal is called.
        system("echo >> stderr.log")
        system("cat stderr.tmp >> stderr.log")
        system("echo  >> stderr.log")
    } else {
        # Run command without capturing stderr
        exitcode <- system2(wcom[1], wcom[2:length(wcom)], env = env)
    }
    # Show error message if command fails
    if (checkexit) {
        if (exitcode != 0) {
            error <- scan("stderr.tmp", what = " ", sep = "\n")
            # Exit code 128+n is typically caused by sending termination signal to process
            if (exitcode > 128) {
                msg <- paste("External program was terminated by the system.\n")
                msg <- paste(msg, "This is usually caused by the system running out of memory.\n\n")
                msg <- paste(msg, "Failed command:\n")
                msg <- paste(msg, command, "\n\n")
                msg <- paste(msg, "Error:\n")
                msg <- paste(msg, paste(error, collapse = "\n"))
                stop(paste("CHIPSTER-NOTE: ", msg))
                # User errors typically cause exit code 1, but this depends on program
            } else {
                msg <- paste("External program did not complete succesfully.\n")
                msg <- paste(msg, "Please check your parameters, input files and input file assignment.\n\n")
                msg <- paste(msg, "Failed command:\n")
                msg <- paste(msg, command, "\n\n")
                msg <- paste(msg, "Error:\n")
                msg <- paste(msg, paste(error, collapse = "\n"))
                system("cat stderr.tmp")
                system("ls -l")
                stop(paste("CHIPSTER-NOTE: ", msg))
            }
        }
    }
}

# Changes the file names in a text file to display names according to chipster-inputs.tsv
#
displayNamesToFile <- function(input.file) {
    input.names <- read_input_names()
    # Escape the characters that are special in a sed pattern or replacement
    pattern <- gsub("([][\\\\.*^$/])", "\\\\\\1", input.names[, 1], perl = TRUE)
    replacement <- gsub("([\\\\&/])", "\\\\\\1", input.names[, 2], perl = TRUE)
    # The script goes through a file, so that the shell doesn't see the names, which can contain
    # spaces and parentheses
    sed.file <- tempfile(fileext = ".sed")
    on.exit(unlink(sed.file))
    writeLines(paste("s/", pattern, "/", replacement, "/", sep = ""), sed.file, useBytes = TRUE)
    runExternal(paste("sed -i -f", sed.file, input.file))
}

# Formats and prints out the command to stdout. Input names are substituted with
# display names to aid readability and help spot input assignment errors.
#
documentCommand <- function(command.string) {
    # Substitute input names. As text, not as regular expressions: names can contain e.g. . ( ) +
    input.names <- tryCatch(read_input_names(), error = function(e) {
        print("no inputs")
        NULL
    })

    if (!is.null(input.names)) {
        for (i in seq_len(nrow(input.names))) {
            command.string <- gsub(input.names[i, 1], input.names[i, 2], command.string, fixed = TRUE)
        }
    }
    cat("##", "COMMAND:", command.string, "\n")
}
