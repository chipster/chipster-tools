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
    report(isTRUE(all.equal(current, target)), info, c(
      paste("  current:", paste(deparse(current), collapse = "")),
      paste("  target: ", paste(deparse(target), collapse = ""))
    ))
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

ok.names <- c("sample_1.fq", "S1+x-2.fastq.gz", utf8_bytes(c(0x6e, 0xc3, 0xa4, 0x2e, 0x66, 0x71)),
  utf8_bytes(c(0x61, 0xcc, 0x88, 0x2e, 0x66, 0x71)))
expect_equal(safe_name(ok.names), ok.names,
  info = "safe_name accepts letters, also non-ASCII and decomposed ones, numbers and . _ + -"
)
expect_equal(safe_name(character(0)), character(0), info = "safe_name accepts an empty vector")
for (bad in c("a b.fq", "a;b.fq", "$(id).fq", "-x.fq", "..", "")) {
  expect_error(safe_name(c("ok.fq", bad)), "CHIPSTER-NOTE: File and sample names",
    info = paste0("safe_name rejects \"", bad, "\"")
  )
}
# The note contains the invalid name, so match it with useBytes: in a UTF-8 locale grepl() can't match
# the message otherwise
msg <- tryCatch({
  safe_name(rawToChar(as.raw(c(0x6e, 0xe4, 0x2e, 0x66, 0x71))))
  "no error"
}, error = function(e) conditionMessage(e))
expect_equal(grepl("CHIPSTER-NOTE: File and sample names", msg, useBytes = TRUE), TRUE,
  info = "safe_name rejects a name that isn't valid UTF-8"
)


setwd(old.wd)
if (run.standalone) {
  cat(if (failures == 0) "all tests passed" else paste(failures, "test(s) FAILED"), "\n")
  if (failures > 0) {
    quit(save = "no", status = 1)
  }
}
