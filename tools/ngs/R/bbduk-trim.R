# TOOL bbduk-trim.R: "Trim QuantSeq reads using BBDuk" (Given a FASTQ file containing Lexogen QuantSeq reads, removes adapters, polyA read-through, and low quality tails. This tool is based on BBDuk.)
# INPUT reads.fq: "FASTQ file" TYPE GENERIC 
# OUTPUT trimmed.fq.gz
# PARAMETER OPTIONAL k: "Kmer length for detecting adapters and polyA" TYPE INTEGER FROM 1 DEFAULT 13 (Contaminants shorter than k will not be found, k must be at least 1.)
# PARAMETER OPTIONAL ktrim: "Trimming direction from a kmer match" TYPE [f: "Don't trim", r: "Trim to the right", l: "Trim to the left"] DEFAULT r (Once a polyA or Truseq adapter kmer is matched in a read, that kmer and all the bases to this direction will be trimmed.)
# PARAMETER OPTIONAL mink: "Minimum length of kmers to report at read tips" TYPE INTEGER DEFAULT 5 (Look for polyA or Truseq adapter kmers down to this length at read tips. 0 means disabled.)
# PARAMETER OPTIONAL qtrim: "Trim read ends based on quality" TYPE [rl: (Trim both ends), f: (Neither end), r: (Right end only), l: (Left end only), w: (Sliding window)] DEFAULT rl (After looking for kmers, remove low quality bases from read ends. Set the quality threshold with the next parameter.)
# PARAMETER OPTIONAL trimq: "Threshold for quality trimming" TYPE DECIMAL DEFAULT 10 (Regions with average base quality below this will be trimmed, if quality trimming is selected. Can be a floating-point number like 7.3.)
# PARAMETER OPTIONAL minlength: "Minimum length for reads to be kept after trimming" TYPE INTEGER DEFAULT 20 (Reads shorter than this after trimming will be discarded.)

source(file.path(chipster.common.path, "tool-utils.R"))
source(file.path(chipster.common.path, "zip-utils.R"))

# Uncompress input if compressed
unzipIfGZipFile("reads.fq")

# File paths
java.path <- paste("PATH=", c(file.path(chipster.tools.path,"rtg","jre","bin")),":$PATH",sep="",collapse="")
bbduk.binary <- c(file.path(chipster.tools.path, "bbmap", "bbduk.sh"))
polya.file <- c(file.path(chipster.tools.path, "bbmap", "resources", "polyA.fa.gz"))
trueseq.file <- c(file.path(chipster.tools.path, "bbmap", "resources", "truseq.fa.gz"))

# Run BBduk
system(paste(bbduk.binary,"-v 2> version.tmp"))
version <- system("grep Version version.tmp",intern = TRUE)
documentVersion("BBDuk",version)
bbduk.options <- paste("-Xmx154m in=reads.fq ", "out=trimmed.fq ", "ref=", polya.file,",",trueseq.file, " k=",k," ktrim=",ktrim," mink=",mink, " qtrim=",qtrim, " trimq=",trimq," minlength=",minlength,sep="",collapse="")
bbduk.command <- paste(bbduk.binary, bbduk.options)
documentCommand(bbduk.command)
#system(bbduk.command)
runExternal(bbduk.command, java.path)

# Compress output
system("gzip trimmed.fq")

# Output names
inputnames <- read_input_definitions()
basename <- strip_name(inputnames$reads.fq)

# Make a matrix of output names
outputnames <- matrix(NA, nrow=1, ncol=2)
outputnames[1,] <- c("trimmed.fq.gz", paste(basename, "_trimmed.fq.gz", sep=""))

# Write output definitions file
write_output_definitions(outputnames)