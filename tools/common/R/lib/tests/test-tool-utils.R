# Tests for the input name helpers in tool-utils.R: read_input_names, make_input_list,
# displayNamesToFile and documentCommand.
#
# Written in the style of the tinytest package, but runs without it on any R in tools-bin. Run it in
# a comp container, e.g. on this VM:
#
#   tool-dev/tool-dev.py shell bowtie2-paired-end.R -- /opt/chipster/tools/R-4.1.1/bin/R --vanilla --slave \
#     -f /opt/chipster/toolbox/tools/common/R/lib/tests/test-tool-utils.R
#
# The R process exits with status 1 if a test fails. Set CHIPSTER_TEST_LARGE=1 to also run the test with a
# 220 MB file.

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
  # all.equal, like tinytest
  expect_equal <- function(current, target, info = "") {
    report(isTRUE(all.equal(current, target)), info, c(
      paste("  current:", paste(deparse(current), collapse = "")),
      paste("  target: ", paste(deparse(target), collapse = ""))
    ))
  }
  expect_true <- function(current, info = "") {
    expect_equal(current, TRUE, info)
  }
  expect_error <- function(current, pattern = ".*", info = "") {
    msg <- tryCatch({
      current
      NULL
    }, error = function(e) conditionMessage(e))
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

# Writes chipster-inputs.tsv like comp does (ToolUtils.writeInputDescription)
write_inputs <- function(dataset.names, input.names = sprintf("reads%03d.fq", seq_along(dataset.names))) {
  writeLines(c(
    "# Chipster dataset description file", "# ", "# INPUT_NAME\tDATASET_NAME",
    paste(input.names, dataset.names, sep = "\t")
  ), "chipster-inputs.tsv", useBytes = TRUE)
}

# Runs displayNamesToFile on the given lines and returns the result
display_names <- function(lines) {
  writeLines(lines, "test.txt", useBytes = TRUE)
  displayNamesToFile("test.txt")
  readLines("test.txt")
}

# Returns what documentCommand prints, without the "## COMMAND:" prefix
document_command <- function(command) {
  sub("^## COMMAND: (.*) $", "\\1", capture.output(documentCommand(command)))
}


## read_input_names

write_inputs(c("001", "010", "NA"))
expect_equal(read_input_names()[, 2], c("001", "010", "NA"), info = "read_input_names keeps names like numbers and NA")


## displayNamesToFile

# Dataset names can contain letters, numbers, space and + - _ : . , ( ) (comp's NAME_PATTERN). The other
# characters here are escaped too, in case the pattern is extended.
write_inputs(c("S1 (1)_1.fq", "näyte+x.fq", "a&b/c\\d [x].fq", "two  spaces.fq"))
expect_equal(
  display_names(c("cmd -1 reads001.fq,reads002.fq -2 reads003.fq reads004.fq reads001.fq", "not me: reads001Xfq")),
  c("cmd -1 S1 (1)_1.fq,näyte+x.fq -2 a&b/c\\d [x].fq two  spaces.fq S1 (1)_1.fq", "not me: reads001Xfq"),
  info = "displayNamesToFile handles special characters and repeated names, and doesn't treat . as a regex"
)

writeBin(c(charToRaw("bad byte: "), as.raw(0xff), charToRaw(" reads001.fq\n")), "test.txt")
displayNamesToFile("test.txt")
expect_equal(
  readBin("test.txt", "raw", 100),
  c(charToRaw("bad byte: "), as.raw(0xff), charToRaw(" S1 (1)_1.fq\n")),
  info = "displayNamesToFile keeps bytes that aren't valid UTF-8"
)

write_inputs(c("001", "010", "NA"))
expect_equal(display_names("-1 reads001.fq -2 reads002.fq -3 reads003.fq"), "-1 001 -2 010 -3 NA",
  info = "displayNamesToFile handles names like numbers and NA"
)

# Display names that are the same as other inputs' names, assigned crosswise
write_inputs(c("reads002.fq", "reads001.fq"))
expect_equal(display_names("-1 reads001.fq -2 reads002.fq"), "-1 reads002.fq -2 reads001.fq",
  info = "displayNamesToFile doesn't change a display name again"
)

# An input name inside another one, as in mothur-trimseqs-uniqueseqs.R
write_inputs(c("sample.fastq", "my.oligos"), c("reads", "reads.oligos"))
expect_equal(
  display_names("trim.seqs(fasta=reads, oligos=reads.oligos) reads reads.oligos"),
  "trim.seqs(fasta=sample.fastq, oligos=my.oligos) sample.fastq my.oligos",
  info = "displayNamesToFile doesn't change an input name inside another one"
)

write_inputs(c("A.txt", "B.txt"), c("other.txt", "input"))
expect_equal(display_names("other.txt input"), "A.txt B.txt",
  info = "displayNamesToFile works with an input name that could be part of a placeholder"
)

# An input name is changed only where it is a file name of its own
write_inputs(c("S1 (1).fq", "my genome"), c("reads001.fq", "reference"))
expect_equal(
  display_names(c(
    "-1 reads001.fq,reads001.fq -x /path/reference (reads001.fq)",
    "reference.fasta reads001.fq.gz x.reads001.fq my_reference reference-1"
  )),
  c(
    "-1 S1 (1).fq,S1 (1).fq -x /path/my genome (S1 (1).fq)",
    "reference.fasta reads001.fq.gz x.reads001.fq my_reference reference-1"
  ),
  info = "displayNamesToFile changes a name next to a separator, but not inside a longer file name"
)
expect_equal(display_names("reads001.fq"), "S1 (1).fq",
  info = "displayNamesToFile changes a name that is the whole line"
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
  writeLines(c(
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\treads001.fq\treads002.fq",
    paste0("chr1\t", 1:n, "\t.\tA\tG\t50\tPASS\tDP=10\tGT\t0/1\t1/1")
  ), "large.vcf")
  seconds <- system.time(displayNamesToFile("large.vcf"))[["elapsed"]]
  cat("     ", round(file.size("large.vcf") / 1e6), "MB in", seconds, "s\n")
  expect_equal(
    readLines("large.vcf", n = 1),
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1 (1).bam\tS2.bam",
    info = "displayNamesToFile handles a large file"
  )
  expect_equal(length(readLines("large.vcf")), n + 1L, info = "displayNamesToFile keeps all lines of a large file")
}


## documentCommand

write_inputs(c("reads002.fq", "reads001.fq"))
expect_equal(
  document_command("bowtie2 -1 reads001.fq -2 reads002.fq reads001Xfq"),
  "bowtie2 -1 reads002.fq -2 reads001.fq reads001Xfq",
  info = "documentCommand doesn't change a display name again, and doesn't treat . as a regex"
)

write_inputs(c("sample.fastq", "my.oligos"), c("reads", "reads.oligos"))
expect_equal(
  document_command("trim.seqs(fasta=reads, oligos=reads.oligos)"),
  "trim.seqs(fasta=sample.fastq, oligos=my.oligos)",
  info = "documentCommand doesn't change an input name inside another one"
)

write_inputs(c("A.txt", "B.txt"), c("other.txt", "input"))
expect_equal(document_command("other.txt input"), "A.txt B.txt",
  info = "documentCommand works with an input name that could be part of a placeholder"
)

write_inputs(c("S1 (1).fq", "my genome"), c("reads001.fq", "reference"))
expect_equal(
  document_command("-1 reads001.fq,reads001.fq -x /path/reference -R reference.fasta reads001.fq.gz"),
  "-1 S1 (1).fq,S1 (1).fq -x /path/my genome -R reference.fasta reads001.fq.gz",
  info = "documentCommand changes a name next to a separator, but not inside a longer file name"
)


## make_input_list

write_inputs(c("S1 (1)_1.fq", "S2_1.fq"))
writeLines(c("S2_1.fq", "S1 (1)_1.fq"), "list.txt")
expect_equal(make_input_list("list.txt"), c("reads002.fq", "reads001.fq"),
  info = "make_input_list finds names with regex characters, in the order of the list"
)

writeBin(charToRaw("S2_1.fq\r\nS1 (1)_1.fq\r\n"), "list.txt")
expect_equal(make_input_list("list.txt"), c("reads002.fq", "reads001.fq"),
  info = "make_input_list reads a list file with Windows line endings"
)

write_inputs(c("S1 ", "S1"))
writeLines("S1 ", "list.txt")
expect_equal(make_input_list("list.txt"), "reads001.fq",
  info = "make_input_list finds a name that ends with a space"
)

write_inputs(c("001", "010", "NA"))
writeLines(c("010", "NA"), "list.txt")
expect_equal(make_input_list("list.txt"), c("reads002.fq", "reads003.fq"),
  info = "make_input_list finds names like numbers and NA"
)

write_inputs(c("S1 (1)_1.fq", "S1 (1)_1.fq", "S2_1.fq"))
writeLines("S2_1.fq", "list.txt")
expect_equal(make_input_list("list.txt"), "reads003.fq",
  info = "make_input_list allows duplicate dataset names that aren't listed"
)
writeLines(c("S2_1.fq", "S1 (1)_1.fq"), "list.txt")
expect_error(make_input_list("list.txt"), "Several selected files have the same name: S1 \\(1\\)_1.fq",
  info = "make_input_list stops when a listed name matches several datasets"
)

write_inputs(c("S1.fq", "S2.fq"))
writeLines(c("S1.fq", "S3.fq"), "list.txt")
expect_error(make_input_list("list.txt"), "has not been selected: S3.fq",
  info = "make_input_list stops when a listed file isn't selected"
)


setwd(old.wd)
if (run.standalone) {
  cat(if (failures == 0) "all tests passed" else paste(failures, "test(s) FAILED"), "\n")
  if (failures > 0) {
    quit(save = "no", status = 1)
  }
}
