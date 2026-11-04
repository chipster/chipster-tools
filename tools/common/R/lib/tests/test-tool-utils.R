# Tests for the helpers in tool-utils.R.
#
# Written in the style of the tinytest package, but runs without it on any R in tools-bin. Run it in
# a comp container with tools-bin mounted (see chipster-openshift/k3s/tool-script-dev.md), on the old
# and the new R:
#
#   /opt/chipster/tools/R-4.1.1/bin/R --vanilla --slave -f /opt/chipster/toolbox/tools/common/R/lib/tests/test-tool-utils.R
#
# The R process exits with status 1 if a test fails.

# Without tinytest, define the expectations used here. Each prints its result and counts the failures.
if (!exists("expect_equal")) {
    failures <- 0
    report <- function(ok, info, details = NULL) {
        cat(if (ok) "ok  " else "FAIL", info, "\n")
        if (!ok) {
            failures <<- failures + 1
            cat(details, sep = "\n")
        }
        invisible(ok)
    }
    # all.equal, like tinytest. An error in the tested code is a failure, not the end of the test run.
    expect_equal <- function(current, target, info = "") {
        current <- tryCatch(current, error = function(e) paste("error:", conditionMessage(e)))
        report(
            isTRUE(all.equal(current, target)),
            info,
            c(
                paste("  current:", paste(deparse(current), collapse = "")),
                paste("  target: ", paste(deparse(target), collapse = ""))
            )
        )
    }
    expect_error <- function(current, pattern = ".*", info = "") {
        msg <- tryCatch(
            {
                current
                NULL
            },
            error = function(e) conditionMessage(e)
        )
        report(!is.null(msg) && grepl(pattern, msg), info, paste("  error:", if (is.null(msg)) "none" else msg))
    }
    run.standalone <- TRUE
} else {
    run.standalone <- FALSE
}

# tool-utils.R is in the parent directory of this file. The file is given as "-f FILE" to R and as
# "--file=FILE" to Rscript. tinytest runs a test file in its own directory.
args <- commandArgs()
file.arg <- c(sub("^--file=", "", grep("^--file=", args, value = TRUE)), args[which(args == "-f") + 1])
tests.dir <- if (length(file.arg) > 0) dirname(file.arg[1]) else getwd()
source(file.path(normalizePath(tests.dir), "..", "tool-utils.R"))

cat("R", paste(R.version$major, R.version$minor, sep = "."), "locale", Sys.getlocale("LC_CTYPE"), "\n")
test.dir <- tempfile("test-tool-utils-")
dir.create(test.dir)
old.wd <- setwd(test.dir)


# --- safe_name --------------------------------------------------------------------------------------

# Builds a string from bytes, marked as UTF-8 like a name from a list file after safe_name's marking
utf8_bytes <- function(bytes) {
    s <- rawToChar(as.raw(bytes))
    Encoding(s) <- "UTF-8"
    s
}

# On R 3.2.3 this tests the fallback of valid_utf8(), on later R validUTF8() itself
expect_equal(
    valid_utf8(c(
        "sample_1.fq",
        utf8_bytes(c(0x6e, 0xc3, 0xa4, 0x2e, 0x66, 0x71)), # nä.fq
        utf8_bytes(c(0x61, 0xcc, 0x88)), # a + combining diaeresis
        utf8_bytes(c(0xf0, 0x9f, 0x98, 0x80)), # 4-byte character
        utf8_bytes(c(0x6e, 0xe4, 0x2e, 0x66, 0x71)), # nä.fq in Latin-1
        utf8_bytes(c(0xc0, 0xaf)), # overlong /
        utf8_bytes(c(0xed, 0xa0, 0x80)), # surrogate
        utf8_bytes(c(0x61, 0xe2, 0x82)), # truncated
        utf8_bytes(c(0xf4, 0x90, 0x80, 0x80)), # above U+10FFFF
        utf8_bytes(c(0x80, 0x61)) # lone continuation byte
    )),
    c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
    info = "valid_utf8 tells valid UTF-8 from invalid"
)
expect_equal(valid_utf8(character(0)), logical(0), info = "valid_utf8 accepts an empty vector")
expect_equal(
    valid_utf8(c(
        utf8_bytes(c(0xef, 0xb7, 0x90)), # U+FDD0, a noncharacter that utf8ToInt() accepts in R 3.2.3
        utf8_bytes(c(0xf4, 0x8f, 0xbf, 0xbf)), # U+10FFFF, the last code point
        utf8_bytes(c(0xe0, 0x80, 0xaf)), # 3-byte overlong /
        utf8_bytes(c(0xf0, 0x80, 0x80, 0xaf)) # 4-byte overlong /
    )),
    c(TRUE, TRUE, FALSE, FALSE),
    info = "valid_utf8 accepts the last code point and rejects longer overlong forms"
)
expect_equal(valid_utf8(NA_character_), TRUE, info = "valid_utf8 treats NA as valid, like validUTF8")

# The known difference of the R 3.2.3 fallback: it rejects the noncharacters U+FFFE and U+FFFF
nonchars <- c(utf8_bytes(c(0xef, 0xbf, 0xbe)), utf8_bytes(c(0xef, 0xbf, 0xbf)))
if (exists("validUTF8", mode = "function", envir = baseenv())) {
    expect_equal(valid_utf8(nonchars), c(TRUE, TRUE), info = "valid_utf8 accepts U+FFFE and U+FFFF (validUTF8)")
} else {
    expect_equal(
        valid_utf8(nonchars),
        c(FALSE, FALSE),
        info = "valid_utf8 rejects U+FFFE and U+FFFF (R 3.2.3 fallback)"
    )
}

ok.names <- c(
    "sample_1.fq",
    "S1+x-2.fastq.gz",
    utf8_bytes(c(0x6e, 0xc3, 0xa4, 0x2e, 0x66, 0x71)),
    utf8_bytes(c(0x61, 0xcc, 0x88, 0x2e, 0x66, 0x71))
)
expect_equal(
    safe_name(ok.names),
    ok.names,
    info = "safe_name accepts letters, also non-ASCII and decomposed ones, numbers and . _ + -"
)
expect_equal(safe_name(character(0)), character(0), info = "safe_name accepts an empty vector")
for (bad in c("a b.fq", "a;b.fq", "$(id).fq", "-x.fq", "..", "")) {
    expect_error(
        safe_name(c("ok.fq", bad)),
        "CHIPSTER-NOTE: File and sample names",
        info = paste0("safe_name rejects \"", bad, "\"")
    )
}
# The note contains the invalid name, so match it with useBytes: in a UTF-8 locale grepl() can't match
# the message otherwise
msg <- tryCatch(
    {
        safe_name(rawToChar(as.raw(c(0x6e, 0xe4, 0x2e, 0x66, 0x71))))
        "no error"
    },
    error = function(e) conditionMessage(e)
)
expect_equal(
    grepl("CHIPSTER-NOTE: File and sample names", msg, useBytes = TRUE),
    TRUE,
    info = "safe_name rejects a name that isn't valid UTF-8"
)

# --- read_list_file ---------------------------------------------------------------------------------

# Writes the given bytes to a file and returns its name
write_bytes <- function(bytes, name = "list.txt") {
    writeBin(bytes, name)
    name
}
bom <- as.raw(c(0xef, 0xbb, 0xbf))

expect_equal(
    read_list_file(write_bytes(c(bom, charToRaw("S2_1.fq\r\n\r\nS2_2.fq\r\n")))),
    c("S2_1.fq", "S2_2.fq"),
    info = "read_list_file removes the BOM, accepts CRLF and skips a blank line"
)
expect_equal(
    read_list_file(write_bytes(charToRaw("a.fq\nb.fq\n\n\n"))),
    c("a.fq", "b.fq"),
    info = "read_list_file skips blank lines at the end"
)
expect_equal(
    read_list_file(write_bytes(charToRaw("a.fq\n  \n\t\n \r\n\f\nb.fq\n"))),
    c("a.fq", "b.fq"),
    info = "read_list_file skips lines with only whitespace"
)
nbsp <- as.raw(c(0xc2, 0xa0))
expect_equal(
    read_list_file(write_bytes(c(charToRaw("a.fq\n"), nbsp, charToRaw("\n "), nbsp, charToRaw("\nb.fq\n")))),
    c("a.fq", "b.fq"),
    info = "read_list_file skips lines with only non-breaking spaces"
)
with.nbsp <- read_list_file(write_bytes(c(charToRaw("a"), nbsp, charToRaw("b.fq\n"))))
expect_equal(
    charToRaw(with.nbsp[1]),
    c(charToRaw("a"), nbsp, charToRaw("b.fq")),
    info = "read_list_file keeps a non-breaking space inside a name"
)
expect_equal(
    read_list_file(write_bytes(c(bom, charToRaw("\na.fq\n")))),
    "a.fq",
    info = "read_list_file skips a first line with only the BOM"
)
expect_equal(
    read_list_file(write_bytes(charToRaw(" a b.fq \nc.fq\n"))),
    c(" a b.fq ", "c.fq"),
    info = "read_list_file keeps spaces in and around names"
)
expect_equal(
    read_list_file(write_bytes(charToRaw("s1\ts1_1.fq\ts1_2.fq\n\ns2\ts2_1.fq\ts2_2.fq\n"))),
    c("s1\ts1_1.fq\ts1_2.fq", "s2\ts2_1.fq\ts2_2.fq"),
    info = "read_list_file keeps tab-separated rows"
)
expect_equal(
    read_list_file(write_bytes(raw(0))),
    character(0),
    info = "read_list_file returns character(0) for an empty file"
)
expect_equal(
    read_list_file(write_bytes(charToRaw("\n\n \n"))),
    character(0),
    info = "read_list_file returns character(0) for a file with only blank lines"
)
# A line that isn't valid UTF-8 (ä in Latin-1) must not stop the function, and comes back as it is
latin1.name <- as.raw(c(0x6e, 0xe4, 0x2e, 0x66, 0x71))
latin1 <- read_list_file(write_bytes(c(latin1.name, charToRaw("\n\nb.fq\n"))))
expect_equal(length(latin1), 2L, info = "read_list_file accepts a line that isn't valid UTF-8")
expect_equal(charToRaw(latin1[1]), latin1.name, info = "read_list_file returns the invalid line unchanged")
expect_equal(latin1[2], "b.fq", info = "read_list_file skips the blank line after the invalid line")

# UTF-16, which e.g. Notepad calls "Unicode", can't be read: stop with a note
utf16le <- as.raw(c(0xff, 0xfe, 0x61, 0x00, 0x2e, 0x00, 0x66, 0x00, 0x71, 0x00, 0x0d, 0x00, 0x0a, 0x00))
expect_error(
    read_list_file(write_bytes(utf16le)),
    "UTF-16",
    info = "read_list_file stops with a note for a UTF-16LE file"
)
utf16be <- as.raw(c(0xfe, 0xff, 0x00, 0x61, 0x00, 0x2e, 0x00, 0x66, 0x00, 0x71, 0x00, 0x0a))
expect_error(
    read_list_file(write_bytes(utf16be)),
    "UTF-16",
    info = "read_list_file stops with a note for a UTF-16BE file"
)

# No final newline. In R 3.2.3, readLines(gzfile(filename)) dropped the last line of such a file.
# The function must not warn about it either, so a warning is returned instead of the lines.
read_quietly <- function(filename) {
    tryCatch(read_list_file(filename), warning = function(w) paste("warning:", conditionMessage(w)))
}
expect_equal(
    read_quietly(write_bytes(charToRaw("a.fq\nb.fq"))),
    c("a.fq", "b.fq"),
    info = "read_list_file keeps the last line of a file without a final newline, without a warning"
)
expect_equal(
    read_quietly(write_bytes(c(bom, charToRaw("S2_1.fq\r\n\r\nS2_2.fq")))),
    c("S2_1.fq", "S2_2.fq"),
    info = "read_list_file keeps the last line after a BOM, CRLF and a blank line"
)

# Gzipped list files, with and without a final newline
con <- gzfile("list.txt.gz", "w")
writeLines(c("a.fq", "", "b.fq"), con)
close(con)
expect_equal(read_list_file("list.txt.gz"), c("a.fq", "b.fq"), info = "read_list_file reads a gzipped file")
con <- gzfile("list.txt.gz", "w")
cat("a.fq\nb.fq", file = con)
close(con)
expect_equal(
    read_quietly("list.txt.gz"),
    c("a.fq", "b.fq"),
    info = "read_list_file keeps the last line of a gzipped file without a final newline"
)

# --- read_input_names, make_input_list and displayNamesToFile -------------------------------------

# Set CHIPSTER_TEST_LARGE=1 to also run the test with a 220 MB file.

# Writes chipster-inputs.tsv like comp does (ToolUtils.writeInputDescription)
write_inputs <- function(dataset.names, input.names = sprintf("reads%03d.fq", seq_along(dataset.names))) {
    writeLines(
        c(
            "# Chipster dataset description file",
            "# ",
            "# INPUT_NAME\tDATASET_NAME",
            paste(input.names, dataset.names, sep = "\t")
        ),
        "chipster-inputs.tsv",
        useBytes = TRUE
    )
}

# Runs displayNamesToFile on the given lines and returns the result
display_names <- function(lines) {
    writeLines(lines, "test.txt", useBytes = TRUE)
    displayNamesToFile("test.txt")
    readLines("test.txt")
}


## read_input_names

write_inputs(c("001", "010", "NA"))
expect_equal(read_input_names()[, 2], c("001", "010", "NA"), info = "read_input_names keeps names like numbers and NA")


## displayNamesToFile

# Dataset names can contain letters, numbers, space and + - _ : . , ( ) (comp's NAME_PATTERN). The other
# characters here are escaped too, in case the pattern is extended.
write_inputs(c("S1 (1)_1.fq", "näyte+x.fq", "a&b/c\\d [x].fq", "two  spaces.fq"))
expect_equal(
    display_names(c("cmd -1 reads001.fq,reads002.fq -2 reads003.fq reads004.fq", "not me: reads001Xfq")),
    c("cmd -1 S1 (1)_1.fq,näyte+x.fq -2 a&b/c\\d [x].fq two  spaces.fq", "not me: reads001Xfq"),
    info = "displayNamesToFile handles special characters, and doesn't treat . as a regex"
)

writeBin(c(charToRaw("bad byte: "), as.raw(0xff), charToRaw(" reads001.fq\n")), "test.txt")
displayNamesToFile("test.txt")
expect_equal(
    readBin("test.txt", "raw", 100),
    c(charToRaw("bad byte: "), as.raw(0xff), charToRaw(" S1 (1)_1.fq\n")),
    info = "displayNamesToFile keeps bytes that aren't valid UTF-8"
)

write_inputs(c("001", "010", "NA"))
expect_equal(
    display_names("-1 reads001.fq -2 reads002.fq -3 reads003.fq"),
    "-1 001 -2 010 -3 NA",
    info = "displayNamesToFile handles names like numbers and NA"
)

write_inputs(sprintf("S%d.fq", 1:12))
expect_equal(
    display_names(paste(sprintf("reads%03d.fq", 1:12), collapse = " ")),
    paste(sprintf("S%d.fq", 1:12), collapse = " "),
    info = "displayNamesToFile works with more than 9 inputs"
)

if (Sys.getenv("CHIPSTER_TEST_LARGE") == "1") {
    write_inputs(c("S1 (1).bam", "S2.bam"))
    n <- 5e6
    writeLines(
        c(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\treads001.fq\treads002.fq",
            paste0("chr1\t", 1:n, "\t.\tA\tG\t50\tPASS\tDP=10\tGT\t0/1\t1/1")
        ),
        "large.vcf"
    )
    seconds <- system.time(displayNamesToFile("large.vcf"))[["elapsed"]]
    cat("     ", round(file.size("large.vcf") / 1e6), "MB in", seconds, "s\n")
    expect_equal(
        readLines("large.vcf", n = 1),
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1 (1).bam\tS2.bam",
        info = "displayNamesToFile handles a large file"
    )
    expect_equal(length(readLines("large.vcf")), n + 1, info = "displayNamesToFile keeps all lines of a large file")
}


## make_input_list

write_inputs(c("S1 (1)_1.fq", "S2_1.fq"))
writeLines(c("S2_1.fq", "S1 (1)_1.fq"), "list.txt")
expect_equal(
    make_input_list("list.txt"),
    c("reads002.fq", "reads001.fq"),
    info = "make_input_list finds names with regex characters, in the order of the list"
)

writeBin(charToRaw("S2_1.fq\r\nS1 (1)_1.fq\r\n"), "list.txt")
expect_equal(
    make_input_list("list.txt"),
    c("reads002.fq", "reads001.fq"),
    info = "make_input_list reads a list file with Windows line endings"
)

write_inputs(c("S1 ", "S1"))
writeLines("S1 ", "list.txt")
expect_equal(make_input_list("list.txt"), "reads001.fq", info = "make_input_list finds a name that ends with a space")

write_inputs(c("001", "010", "NA"))
writeLines(c("010", "NA"), "list.txt")
expect_equal(
    make_input_list("list.txt"),
    c("reads002.fq", "reads003.fq"),
    info = "make_input_list finds names like numbers and NA"
)

write_inputs(c("S1 (1)_1.fq", "S1 (1)_1.fq", "S2_1.fq"))
writeLines("S2_1.fq", "list.txt")
expect_equal(
    make_input_list("list.txt"),
    "reads003.fq",
    info = "make_input_list allows duplicate dataset names that aren't listed"
)
writeLines(c("S2_1.fq", "S1 (1)_1.fq"), "list.txt")
expect_error(
    make_input_list("list.txt"),
    "Several selected files have the same name: S1 \\(1\\)_1.fq",
    info = "make_input_list stops when a listed name matches several datasets"
)

write_inputs(c("S1.fq", "S2.fq"))
writeLines(c("S1.fq", "S3.fq"), "list.txt")
expect_error(
    make_input_list("list.txt"),
    "has not been selected: S3.fq",
    info = "make_input_list stops when a listed file isn't selected"
)

# The notes name all the files once, in one sentence
note_of <- function(expr) {
    tryCatch(
        {
            expr
            "no error"
        },
        error = function(e) conditionMessage(e)
    )
}
writeLines(c("S1.fq", "S3.fq", "S4.fq"), "list.txt")
expect_equal(
    note_of(make_input_list("list.txt")),
    "CHIPSTER-NOTE:  Input file list includes one or more files that has not been selected: S3.fq, S4.fq",
    info = "make_input_list names all the files that aren't selected in one note"
)
writeLines(c("S1.fq", "S2.fq", "S1.fq", "S2.fq", "S1.fq"), "list.txt")
expect_equal(
    note_of(make_input_list("list.txt")),
    paste(
        "CHIPSTER-NOTE:  Input file list contains one or more duplicated entries: S1.fq, S2.fq",
        "\nFilenames must be unique so they can be correctly assigned."
    ),
    info = "make_input_list names each duplicated entry once, in one note"
)

writeBin(raw(0), "list.txt")
expect_error(
    make_input_list("list.txt"),
    "CHIPSTER-NOTE: The input file list list.txt is empty",
    info = "make_input_list stops with a note for an empty list"
)
writeLines(c("", "  ", ""), "list.txt")
expect_error(
    make_input_list("list.txt"),
    "CHIPSTER-NOTE: The input file list list.txt is empty",
    info = "make_input_list stops with a note for a list with only blank lines"
)
write_inputs(c("S1.fq", "my reads.txt"), c("reads001.fq", "list.txt"))
expect_error(
    make_input_list("list.txt"),
    "CHIPSTER-NOTE: The input file list my reads.txt is empty",
    info = "make_input_list names an empty list by its dataset name"
)
write_inputs(c("S1.fq", "S2.fq"))

# read_list_file() does the reading: byte order mark, blank lines and a missing final newline
writeBin(c(as.raw(c(0xef, 0xbb, 0xbf)), charToRaw("S2.fq\r\n\r\nS1.fq")), "list.txt")
expect_equal(
    make_input_list("list.txt"),
    c("reads002.fq", "reads001.fq"),
    info = "make_input_list reads a list with a BOM, a blank line and no final newline"
)
writeBin(as.raw(c(0xff, 0xfe, 0x53, 0x00, 0x31, 0x00, 0x2e, 0x00, 0x66, 0x00, 0x71, 0x00, 0x0a, 0x00)), "list.txt")
expect_error(make_input_list("list.txt"), "UTF-16", info = "make_input_list stops with a note for a UTF-16 list")


## documentCommand

command_of <- function(command) {
    out <- capture.output(documentCommand(command))
    out[length(out)]
}
# Dataset names can have only the characters of comp's ToolUtils.NAME_PATTERN
write_inputs(c("001", "S1 (1)+x.fq", "NA"))
expect_equal(
    command_of("bowtie2 -1 reads001.fq -2 reads002.fq -U reads003.fq"),
    "## COMMAND: bowtie2 -1 001 -2 S1 (1)+x.fq -U NA ",
    info = "documentCommand shows names like numbers or NA, and names with regex characters"
)
expect_equal(
    command_of("cat reads001Xfq"),
    "## COMMAND: cat reads001Xfq ",
    info = "documentCommand doesn't treat the . of an input name as any character"
)
unlink("chipster-inputs.tsv")
expect_equal(
    capture.output(documentCommand("ls")),
    c("[1] \"no inputs\"", "## COMMAND: ls "),
    info = "documentCommand prints the command as it is without chipster-inputs.tsv"
)

setwd(old.wd)
if (run.standalone) {
    cat(if (failures == 0) "all tests passed" else paste(failures, "test(s) FAILED"), "\n")
    if (failures > 0) {
        quit(save = "no", status = 1)
    }
}
