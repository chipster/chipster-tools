# Tests for the helpers in vcf-utils.R.
#
# Written in the style of the tinytest package, but runs without it on any R in tools-bin. Run it in
# a comp container with tools-bin mounted (see chipster-openshift/k3s/tool-script-dev.md), on the old
# and the new R:
#
#   /opt/chipster/tools/R-4.1.1/bin/R --vanilla --slave -f /opt/chipster/toolbox/tools/common/R/lib/tests/test-vcf-utils.R
#
# The R process exits with status 1 if a test fails.

# Without tinytest, define the expectations used here. Each prints its result and counts the failures.
if (!exists("expect_equal")) {
    failures <- 0
    # all.equal, like tinytest. An error in the tested code is a failure, not the end of the test run.
    expect_equal <- function(current, target, info = "") {
        current <- tryCatch(current, error = function(e) paste("error:", conditionMessage(e)))
        ok <- isTRUE(all.equal(current, target))
        cat(if (ok) "ok  " else "FAIL", info, "\n")
        if (!ok) {
            failures <<- failures + 1
            cat(paste("  current:", paste(deparse(current), collapse = "")), "\n")
            cat(paste("  target: ", paste(deparse(target), collapse = "")), "\n")
        }
        invisible(ok)
    }
    run.standalone <- TRUE
} else {
    run.standalone <- FALSE
}

# The libraries are in the parent directory of this file. The file is given as "-f FILE" to R and as
# "--file=FILE" to Rscript. tinytest runs a test file in its own directory.
args <- commandArgs()
file.arg <- c(sub("^--file=", "", grep("^--file=", args, value = TRUE)), args[which(args == "-f") + 1])
tests.dir <- if (length(file.arg) > 0) dirname(file.arg[1]) else getwd()
# zip-utils.R sources tool-utils.R from chipster.common.lib.path, which comp sets in a job.
chipster.common.lib.path <- normalizePath(file.path(tests.dir, ".."))
source(file.path(chipster.common.lib.path, "zip-utils.R"))
source(file.path(chipster.common.lib.path, "vcf-utils.R"))

cat("R", paste(R.version$major, R.version$minor, sep = "."), "\n")
test.dir <- tempfile("test-vcf-utils-")
dir.create(test.dir)
old.wd <- setwd(test.dir)

header <- c(
    "##fileformat=VCFv4.2",
    "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">",
    "##contig=<ID=1,length=3000>",
    "##contig=<ID=10,length=3000>",
    "##contig=<ID=X,length=3000>",
    "##contig=<ID=MT,length=3000>",
    "##contig=<ID=GL000192.1,length=3000>",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
)
records <- c(
    "1\t500\t.\tA\tC\t.\t.\tAF=0.1",
    "10\t800\t.\tA\t<DEL>\t.\t.\tSVTYPE=DEL;END=900",
    "X\t700\t.\tA\tG\t.\t.\tHGVS=c.700A>G",
    "MT\t100\t.\tA\tG\t.\t.\t.",
    "GL000192.1\t100\t.\tA\tG\t.\t.\t."
)
writeLines(c(header, records), "nochr.vcf")

addChrToVCF("nochr.vcf", "chr.vcf")
chr.vcf <- readLines("chr.vcf")
expect_equal(
    chr.vcf[grep("^#", chr.vcf, invert = TRUE)],
    c(
        "chr1\t500\t.\tA\tC\t.\t.\tAF=0.1",
        "chr10\t800\t.\tA\t<DEL>\t.\t.\tSVTYPE=DEL;END=900",
        "chrX\t700\t.\tA\tG\t.\t.\tHGVS=c.700A>G",
        "chrMT\t100\t.\tA\tG\t.\t.\t.",
        "GL000192.1\t100\t.\tA\tG\t.\t.\t."
    ),
    info = "addChrToVCF adds chr to the records of numbered chromosomes and X, Y, Z, W, M, and nothing else"
)
expect_equal(
    grep("^##contig", chr.vcf, value = TRUE),
    c(
        "##contig=<ID=chr1,length=3000>",
        "##contig=<ID=chr10,length=3000>",
        "##contig=<ID=chrX,length=3000>",
        "##contig=<ID=chrMT,length=3000>",
        "##contig=<ID=GL000192.1,length=3000>"
    ),
    info = "addChrToVCF adds chr to the IDs of the ##contig header lines"
)
expect_equal(
    chr.vcf[grepl("^#", chr.vcf) & !grepl("^##contig", chr.vcf)],
    header[!grepl("^##contig", header)],
    info = "addChrToVCF leaves the other header lines as they are"
)
expect_equal(getVCFNames("nochr.vcf"), "1", info = "getVCFNames sees no chr before addChrToVCF")
expect_equal(getVCFNames("chr.vcf"), "chr1", info = "getVCFNames sees chr after addChrToVCF")

removeChrFromVCF("chr.vcf", "chr-removed.vcf")
expect_equal(readLines("chr-removed.vcf"), c(header, records), info = "removeChrFromVCF reverses addChrToVCF")

setwd(old.wd)
if (run.standalone) {
    cat(if (failures == 0) "all tests passed" else paste(failures, "test(s) FAILED"), "\n")
    if (failures > 0) {
        quit(save = "no", status = 1)
    }
}
