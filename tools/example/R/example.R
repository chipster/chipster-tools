# TOOL example.R: "R template for tools" (R template for tools)
# INPUT input.tsv: "Required input" TYPE GENERIC (Required input file, any type)
# INPUT OPTIONAL optional.tsv: "Optional input" TYPE GENERIC (Demonstrates OPTIONAL prefix)
# INPUT META phenodata.tsv: "Phenodata" TYPE GENERIC (Demonstrates META prefix; rendered as phenodata in the UI)
# INPUT files{...}.tsv: "Multiple files" TYPE GENERIC (Demonstrates multi-input binding with {...})
# INPUT text.txt: "Plain text" TYPE TEXT
# INPUT strict.txt: "Strict .txt only" TYPE TEXT_STRICT
# INPUT archive.tar: "Tar archive" TYPE TAR
# INPUT reads.bam: "Aligned reads" TYPE BAM
# INPUT sequences.fa: "FASTA sequences" TYPE FASTA
# INPUT reads.fq: "FASTQ reads" TYPE FASTQ
# INPUT annotations.gtf: "Gene annotations" TYPE GTF
# INPUT cdna.tsv: "Two-color microarray" TYPE CDNA
# INPUT affy.cel: "Affymetrix CEL" TYPE AFFY
# INPUT expression.tsv: "Gene expression matrix" TYPE GENE_EXPRS
# INPUT genes.tsv: "Gene list" TYPE GENELIST
# INPUT phenodata-typed.tsv: "Phenodata content type" TYPE PHENODATA
# INPUT oligos.txt: "Mothur oligos" TYPE MOTHUR_OLIGOS
# INPUT names.txt: "Mothur names" TYPE MOTHUR_NAMES
# INPUT groups.txt: "Mothur groups" TYPE MOTHUR_GROUPS
# INPUT stability.files: "Mothur stability" TYPE MOTHUR_STABILITY
# INPUT count.txt: "Mothur count" TYPE MOTHUR_COUNT
# OUTPUT output.tsv
# OUTPUT named.tsv: "Output with label" (Output with a quoted label and description)
# OUTPUT OPTIONAL optional.txt: "Optional output" (Demonstrates OPTIONAL prefix; not always created)
# OUTPUT META phenodata-out.tsv: "Phenodata output" (Demonstrates META prefix; written as phenodata)
# OUTPUT files{...}.txt: "Multiple files" (Demonstrates multi-output binding; create any number of matching files)
# PARAMETER OPTIONAL int: Integer TYPE INTEGER FROM -100 TO 100 DEFAULT 0
# PARAMETER OPTIONAL dec: Decimal TYPE DECIMAL FROM -0.1 TO 1.0 DEFAULT 0
# PARAMETER OPTIONAL string: String TYPE STRING DEFAULT "default"
# PARAMETER OPTIONAL unchecked: "Unchecked string" TYPE UNCHECKED_STRING DEFAULT "default"
# PARAMETER OPTIONAL enum: Enum TYPE [first: First, second: Second, third: Third] DEFAULT first
# PARAMETER OPTIONAL column: Column TYPE COLUMN_SEL
# PARAMETER OPTIONAL metacolumn: "Phenodata column" TYPE METACOLUMN_SEL DEFAULT group
# IMAGE comp-r-4-5-1
# RUNTIME R-4.5.1
# SLOTS 1
# STORAGE 10
# TOOLS_BIN ""


# Extra header directives
# ------------------------------------------------------------------------------
# Beyond TOOL / INPUT / OUTPUT / PARAMETER, the SADL header supports the five
# optional directives below. They must appear in the listed order; the parser
# only accepts them in IMAGE -> RUNTIME -> SLOTS -> STORAGE -> TOOLS_BIN
# sequence and rejects the file with a parse error otherwise. The required
# order is alphabetical, which is the easiest way to remember it.
#
# IMAGE <name>      Override the container image. Use when the tool needs an
#                   image other than the runtime's default (e.g. one with extra
#                   system libraries baked in).
#
# RUNTIME <name>    Select the interpreter and base environment. Real values
#                   include R-4.5.1, R-4.4.3-mothur, python3, etc. Determines
#                   the R/Python version and the set of pre-installed packages
#                   the tool runs against.
#
# SLOTS <int>       Resource request. One slot maps to a per-slot CPU and
#                   memory amount set by the deployment config; SLOTS 5 means
#                   the scheduler reserves 5x the CPU and 5x the memory. The
#                   chipster.threads.max and chipster.memory.max variables
#                   below reflect the resulting limits.
#
# STORAGE <int>     Scratch disk space request (integer, GB). Set this
#                   for tools that produce large intermediate files.
#
# TOOLS_BIN ""      Opt out of mounting the chipster-tools binaries into the
#                   job. Pure R/Python tools that don't shell out to external
#                   binaries set this to "". In practice, the directive is never 
#                   used to point at an alternate path.


# Chipster variables
# ------------------------------------------------------------------------------
# The Chipster runtime injects these variables into the R session before the
# script runs. Use them to construct paths to pre-installed binaries and
# library scripts, or to query the resources allocated to this job. Example
# values shown.

# Path to pre-installed binaries:
chipster.tools.path     # "/opt/chipster-web-server/tools"

# Paths to library script directories:
chipster.common.path     # "../toolbox/tools/common/R"
chipster.common.lib.path # "../toolbox/tools/common/R/lib"
chipster.module.path     # "../toolbox/tools/misc"
chipster.module.lib.path # "../toolbox/tools/misc/lib"
chipster.java.libs.path  # "/opt/chipster-web-server/lib"

# Resources allocated to this job:
chipster.threads.max     # "2"      CPU cores
chipster.memory.max      # "8192"   memory in MB


# Use a library script file
# ------------------------------------------------------------------------------
# In addition to standard R packages (see the next section), the Chipster
# project maintains its own R helper scripts in tools/common/R/lib — e.g.
# tool-utils.R, zip-utils.R, version-utils.R. These are plain R files
# versioned in this repo, not installed packages; they hold helpers written
# by tool authors and shared across many tools. Add new helpers there when
# something is worth reusing.
#
# Load a library script with source(), using chipster.common.lib.path to
# find the directory. Since source() evaluates the file into the current
# environment, the helpers become available with no namespace prefix.
source(file.path(chipster.common.lib.path, "tool-utils.R"))

# Use a function from the library script:
is.fasta <- isFasta("input.tsv")


# Load an R package
# ------------------------------------------------------------------------------
# Most tools need one or more R packages beyond base R. Load them with
# library(). The packages must be pre-installed in the runtime selected by the
# RUNTIME directive — see chipster-defaults.yaml for what's bundled in each.
# Wrap the call in suppressPackageStartupMessages() to keep the job log clean.
suppressPackageStartupMessages(library(ggplot2))


# Validate inputs
# ------------------------------------------------------------------------------
# fileOk() / fileNotOk() / fileCheck() (from tool-utils.R) test that a file
# exists and optionally meets a minimum size or line count. fileCheck() halts
# the job with a CHIPSTER-NOTE on failure; the other two return TRUE / FALSE
# so the script can branch on the result.
fileCheck("input.tsv")                  # required, must exist
fileCheck("input.tsv", minsize = 100)   # must be at least 100 bytes
if (fileNotOk("optional.tsv")) {
    stop(paste("CHIPSTER-NOTE:", "Optional input is empty or missing"))
}


# Notify user with Chipster note
# ------------------------------------------------------------------------------
# Halt the script with a user-facing error when something is wrong with the
# inputs or parameters that the user can fix. The "CHIPSTER-NOTE:" prefix on
# the stop() message makes Chipster render it as a friendly error in the UI
# rather than a stack trace, so the text should be actionable.
if (param.count < 2) {
    stop(paste("CHIPSTER-NOTE:", "Count needs to greater than 2"))
}


# Read phenodata
# ------------------------------------------------------------------------------
# The META phenodata input is a TSV with one row per sample and columns the
# user fills in via the UI (group label, condition, batch, etc.). Read it as
# a data frame and look up the column the user selected with a METACOLUMN_SEL
# parameter (here `metacolumn`, declared in the PARAMETER block above).
phenodata <- read.table("phenodata.tsv", sep = "\t", header = TRUE, check.names = FALSE)
groups <- phenodata[[metacolumn]]


# Iterate over multi-input files
# ------------------------------------------------------------------------------
# A multi-input declaration like `# INPUT files{...}.tsv` matches any number
# of files sharing the declared prefix and suffix. At runtime they appear in
# the job directory with the `{...}` replaced by per-file tokens. Use
# list.files() with a regex to enumerate them; cross-reference chipster-
# inputs.tsv if you need the user-visible display names.
input.files <- list.files(".", pattern = "^files.*\\.tsv$")
for (f in input.files) {
    # process each file
}


# Read input definitions
# ------------------------------------------------------------------------------
# read_input_definitions() (from tool-utils.R) returns a named list mapping
# each declared SADL input name (e.g. "input.tsv") to a vector of values read
# from chipster-inputs.tsv; the first element is the user-visible display name
# of the file the user selected. Use it to log inputs by their original
# filename or to construct output names from the input name.
input.names <- read_input_definitions()
input.display.name <- input.names[["input.tsv"]][[1]]


# Rename outputs
# ------------------------------------------------------------------------------
# By default the user sees the label declared in each OUTPUT directive (or the
# filename if no label was given). To assign display names dynamically, write
# a mapping from the internal output filename to a user-visible dataset name
# into chipster-outputs.tsv. The helper write_output_definitions() (from
# tool-utils.R) writes this file; each row is <internal_name>\t<display_name>.
#
# The snippet below renames outputs declared by these header directives:
#   # OUTPUT output.tsv
#   # OUTPUT files{...}.txt
outputnames <- matrix(NA, nrow = 2, ncol = 2)
outputnames[1, ] <- c("output.tsv", "result-table.tsv")
outputnames[2, ] <- c("files0.txt", "first-extra.txt")
write_output_definitions(outputnames)


# Run a command-line tool
# ------------------------------------------------------------------------------
# Most R tool scripts are thin wrappers around an external binary. Build the
# command from these pieces:
#   - chipster.tools.path:  location of pre-installed binaries
#   - chipster.threads.max: number of CPU cores allocated to this job; pass it
#                           to multi-threaded tools (-t / -p / --threads / -@)
#   - runExternal():        helper from tool-utils.R that captures stderr to
#                           stderr.log, checks the exit code, and aborts the
#                           job with a CHIPSTER-NOTE on failure (preferred over
#                           a raw system() call)
#   - documentCommand():    optional helper that prints the command with input
#                           names substituted for user-visible display names,
#                           useful for the job log
samtools.binary <- file.path(chipster.tools.path, "samtools", "bin", "samtools")
command <- paste(samtools.binary, "sort", "-@", chipster.threads.max,
                 "-o", "output.bam", "input.bam")
documentCommand(command)
runExternal(command)


# Document application version
# ------------------------------------------------------------------------------
# Record the version of any external tool the script invokes. Versions are
# stored alongside the job in the session db and shown in the job details,
# making results reproducible.
#
# documentVersion() comes from version-utils.R, which is sourced automatically
# at the start of every R job — no need to source it. The version string can
# span multiple lines.
version <- system(paste(samtools.binary, "--version | head -1 | cut -d ' ' -f 2"), intern = TRUE)
documentVersion("samtools", version)


# Decompress a gzipped input
# ------------------------------------------------------------------------------
# Input files may arrive gzipped. The zip-utils.R library provides
# unzipIfGZipFile() which detects gzip compression with the unix `file` command
# and decompresses the file in place, leaving the original filename intact (so
# the rest of the script can use the declared input name without change).
source(file.path(chipster.common.lib.path, "zip-utils.R"))
unzipIfGZipFile("input.tsv")


